! *****************************************************************************/
! * License
! *        This file is part of EOF-Library.
! *
! *        EOF-Library is free software: you can redistribute it and/or modify it
! *        under the terms of the GNU General Public License as published by
! *        the Free Software Foundation, either version 3 of the License, or
! *        (at your option) any later version.
! *
! *        EOF-Library is distributed in the hope that it will be useful, but WITHOUT
! *        ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! *        FITNESS FOR A PARTICULAR PURPOSE.    See the GNU General Public License
! *        for more details.
! *
! *        You should have received a copy of the GNU General Public License
! *        along with EOF-Library.    If not, see <http://www.gnu.org/licenses/>.
! *
! ******************************************************************************
! *
! *    Authors: Juris Vencels (EOF Consulting, Latvia)
! *                     Peter Råback (CSC - IT Center for Science, Finland)
! *    
! *    Web:         https://eof-library.com
! *
! *    Original Date: 29.09.2016
! *
! *****************************************************************************/

! procedure see:
! J. Vencels, P. Råback, und V. Geža, „EOF-Library: Open-source Elmer FEM and OpenFOAM 
! coupler for electromagnetics and fluid dynamics“, SoftwareX, Bd. 9, S. 68–72, 
! Jan. 2019, doi: 10.1016/j.softx.2019.01.007.

!------------------------------------------------------------------------------
!> Module for sending interpolated fields to OpenFOAM.
!------------------------------------------------------------------------------
MODULE Elmer2OpenFOAMSolverUtils

    USE DefUtils
    USE MeshUtils
    USE ElementUtils
    USE ParticleUtils

    IMPLICIT NONE

    TYPE OFproc_t
        INTEGER :: reqSend, reqRecv, globalRank
        TYPE(Variable_t), POINTER :: OFVar
        TYPE(Mesh_t), POINTER :: OFMesh
        LOGICAL,POINTER :: foundCells(:)
        INTEGER,POINTER :: foundCellsIndx(:)
        INTEGER :: nFoundCells
        LOGICAL :: boxOverlap
    END TYPE OFproc_t

    ! IS BODY FOR MULTIPLE OF BODIES TO ONE ELMER BODY????


    TYPE(OFproc_t), ALLOCATABLE, TARGET :: OFp(:,:) ! [OFrank][body]
    INTEGER :: totOFRanks, OFRanksStart, ElmerRanksStart, myGlobalRank, &
                                     totGlobalRanks, myLocalRank, totLocalRanks, nVars, &
                                     nBodiesToComm
    REAL(KIND=dp) :: myBoundBox(3,2) ! [x,y,z][min,max]
    REAL(KIND=dp), POINTER :: ELboundBoxes(:,:,:,:) ! [x,y,z][min,max][rank][body]
    INTEGER, POINTER :: OF_EL_overlap(:,:,:) ! [ELrank][OFrank][body]

    INTEGER, ALLOCATABLE :: USE_TARGET_VARS_IN_BODY(:,:)
END MODULE Elmer2OpenFOAMSolverUtils

!------------------------------------------------------------------------------
SUBROUTINE MPI_TEST_SLEEP( req, ierr )

    USE ISO_C_BINDING, ONLY : C_LONG
    USE Elmer2OpenFOAMSolverUtils

    IMPLICIT NONE

    INTERFACE
        SUBROUTINE usleep(n) bind(C)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(c_int32_t), VALUE :: n
        END SUBROUTINE usleep
    END INTERFACE

    !------------------------------------------------------------------------------
    INTEGER :: req, ierr
    LOGICAL :: Flag

    DO WHILE ( .TRUE. )
        CALL MPI_TEST( req, Flag, MPI_STATUS_IGNORE, ierr )
        IF (Flag) EXIT
        CALL usleep(1000_c_int32_t)
    END DO

END SUBROUTINE MPI_TEST_SLEEP

!------------------------------------------------------------------------------
SUBROUTINE findOverlappingBoxes(s)

    USE Elmer2OpenFOAMSolverUtils

    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER :: s
    INTEGER :: ierr,    i
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_ALLGATHER(myBoundBox, 6, MPI_DOUBLE, ELboundBoxes(:,:,:,s), 6, MPI_DOUBLE, ELMER_COMM_WORLD, ierr)

    IF ( myLocalRank==0 ) THEN
        CALL MPI_SEND(ELboundBoxes(:,:,:,s), totLocalRanks*2*3, MPI_DOUBLE, &
                                    OFp(0,s) % globalRank, 1001, MPI_COMM_WORLD, ierr)
        CALL MPI_RECV(OF_EL_overlap(:,:,s), totOFRanks*totLocalRanks, MPI_INTEGER, &
                                    OFp(0,s) % globalRank, 1002, MPI_COMM_WORLD, status, ierr)
    END IF

    CALL MPI_Bcast(OF_EL_overlap(:,:,s), totOFRanks*totLocalRanks, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

    DO i=0,totOFRanks-1
        OFp(i,s) % boxOverlap = (OF_EL_overlap(myLocalRank,i,s)==1)
    END DO

END SUBROUTINE findOverlappingBoxes

!------------------------------------------------------------------------------
SUBROUTINE Elmer2OpenFOAMSolver( Model,Solver,dt,TransientSimulation )
    
    USE DefUtils
    USE Interpolation
    USE MeshUtils
    USE ElementUtils
    USE ParticleUtils
    USE Elmer2OpenFOAMSolverUtils

    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt, exitcond
    LOGICAL :: TransientSimulation

    ! local variables
    !------------------------------------------------------------------------------    
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName
    INTEGER :: i, j, ierr, OFstatus, s
    INTEGER :: status(MPI_STATUS_SIZE)
    LOGICAL :: Found
    REAL(KIND=dp) :: commTime
    CHARACTER(LEN=15) :: timeStr
    INTEGER, POINTER :: Blist(:)

    INTEGER, POINTER :: Varlist(:)
    
    INTERFACE
        SUBROUTINE InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, NewVariables, &
                UseQuadrantTree, Projector, MaskName, FoundNodes, NewMaskPerm, KeepUnfoundNodes )
            USE Types
            TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
            TYPE(Mesh_t), TARGET    :: OldMesh, NewMesh
            LOGICAL, OPTIONAL :: UseQuadrantTree,FoundNodes(:)
            CHARACTER(LEN=*),OPTIONAL :: MaskName
            TYPE(Projector_t), POINTER, OPTIONAL :: Projector
            INTEGER, OPTIONAL, POINTER :: NewMaskPerm(:)    !< Mask the new variable set by the given MaskName when trying to define the interpolation.
             LOGICAL, OPTIONAL :: KeepUnfoundNodes    !< Do not disregard unfound nodes from projector
        END SUBROUTINE InterpolateMeshToMeshQ
    END INTERFACE

    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, SAVE :: VISITED = .FALSE. ! variable anlegen

    !------------------------------------------------------------------------------    

    CALL Info('Elmer2OpenFOAMSolver','-----------------------------------------', Level=4 )

    commTime = MPI_WTIME()

    ! The variable containing the field contributions
    !--------------------------------------------------------------------------
    CALL Info('Elmer2OpenFOAMSolver', 'GetSolverParams()', Level=4 )
    Params => GetSolverParams()
    CALL Info('Elmer2OpenFOAMSolver', 'GetMesh()', Level=4 )
    Mesh => GetMesh()

    CALL Info('Elmer2OpenFOAMSolver', TRIM(Mesh % Name) , Level=4 )
    CALL Info('Elmer2OpenFOAMSolver', 'Reading target variable', Level=4 )

    
    CALL Info('Elmer2OpenFOAMSolver', TRIM(I2S(Mesh % Elements(1) % bodyid)) , Level=4 )
    CALL Info('Elmer2OpenFOAMSolver', TRIM(I2S(Mesh % Elements(2000) % bodyid)) , Level=4 )
    ! maybe not goes in here seconds time that it is called maybe visited is global or
    ! saved during runs?
    
    
    ! basically two possbilities define bodies for each variables 
    ! or
    ! make solver callable multiple times -> by solver name ... 
    ! -> most practical?

    ! TODO (Future):
    ! Add Joule heat conservation ...
    ! Add Target bodies to reduce target mesh size ...

    ! maybe with getbody from defutils? 
    ! TODO: Extend for the definition of bodies that should be coupled to openfoam
    IF (.NOT. VISITED) THEN
        nVars = 0

        ! what are nBodiesToComm ?

        IF( ListCheckPresent( Params,'Bodies') ) THEN
            ! body numbers itself seems to be irrelevant -> this is due tu mesh structure? 
            ! nodes are listet in mesh % nodes
            ! but body ids are defined in Mesh % Elements(:) % bodyid
            ! elements have pointer list for nodes but filtering mesh nodes with body id 
            ! is not implemented yet ...
            BList => ListGetIntegerArray( Params, 'Bodies', Found ) 
            nBodiesToComm = SIZE(BList)
        ELSE
            nBodiesToComm = 1
        END IF

        ! Max 20 Target Variables (20 = nexp in SOLVER.KEYWORDS ...)
        DO i=1,20
            VarName = ListGetString( Params, 'Target Variable '//TRIM(I2S(i)), Found )
            IF(.NOT. Found ) THEN
                IF (i==1) CALL Fatal('Elmer2OpenFOAMSolver','> Target Variable 1 < must exist for the solver!')
                EXIT
            ELSE
             ! Test that the variable exists in the primary mesh
                Var => VariableGet(Mesh % Variables, VarName )
                IF(.NOT. ASSOCIATED( Var ) ) THEN
                    CALL Fatal('Elmer2OpenFOAMSolver','Variable '//TRIM(VarName)//' does not exist in Elmer mesh!')
                ELSE
                    nVars = nVars + 1
                END IF
            END IF
        END DO

        IF(nBodiesToComm > 1) THEN
            ALLOCATE(USE_TARGET_VARS_IN_BODY(nVars, nBodiesToComm)) ! frotran has variable indexing standard is 1-n
            DO s = 1, nBodiesToComm
                ! If not defined assume all variables are target variables for each body ...
                ! zeros equal to do not use that target variable in that body ... 
                Varlist => ListGetIntegerArray( Params, 'Body '//TRIM(I2S(s))//' Use Target Variable', Found )

                IF(.NOT. Found ) THEN
                    ! check indexing ????
                    DO i = 1, nVars
                        USE_TARGET_VARS_IN_BODY(i,s) = 1
                    END DO
                ELSE
                    CALL Info('Elmer2OpenFOAMSolver','Read Body '//TRIM(I2S(s))//' Use Target Variable',Level=4)
                    IF(SIZE(Varlist) == nVars) THEN
                        USE_TARGET_VARS_IN_BODY(:,s) = Varlist
                    ELSE
                        ! FATAL 
                        CALL Fatal('Elmer2OpenFOAMSolver','Body'//TRIM(I2S(s))//'Target Variables List error!')
                    END IF

                END IF
                CALL Info('Elmer2OpenFOAMSolver','Target Variables for Body '//TRIM(I2S(s))//':',Level=4)
                DO i = 1, nVars
                    CALL Info('Elmer2OpenFOAMSolver',TRIM(I2S(USE_TARGET_VARS_IN_BODY(i,s))),Level=4)
                END DO
            END DO
        ELSE
            ! check indexing ????
            DO i = 1, nVars
                USE_TARGET_VARS_IN_BODY(i,s) = 1
            END DO
        END IF


        CALL Info('Elmer2OpenFOAMSolver','Number of target variables: '//TRIM(I2S(nVars)),Level=3)

        ! MPI coupling
        !------------------------------------------------------------------------
        myLocalRank     = ParEnv % MyPE
        totLocalRanks = ParEnv % PEs

        CALL MPI_COMM_RANK( MPI_COMM_WORLD, myGlobalRank, ierr )
        CALL MPI_COMM_SIZE( MPI_COMM_WORLD, totGlobalRanks, ierr )

        totOFRanks = totGlobalRanks - totLocalRanks

        IF(totOFRanks==0) THEN
            CALL Fatal('Elmer2OpenFOAMSolver','MPI communicator does not have OpenFOAM procs!')
        END IF

        IF (myLocalRank == 0) ElmerRanksStart = myGlobalRank
        CALL MPI_BCAST(ElmerRanksStart, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

        IF (ElmerRanksStart==0) THEN
             OFRanksStart = totLocalRanks
        ELSE
             OFRanksStart = 0
        END IF

        CALL Info('Elmer2OpenFOAMSolver', 'Allocating OpenFOAM data structures',Level=3)

        ! nBodiesToComm -> seems to be refered to openfoam bodies? ...
        ALLOCATE( OFp(0:totOFRanks-1, nBodiesToComm) )

        ! Get getboundBox
        myBoundBox(1,1) = MINVAL(Mesh % Nodes % x)
        myBoundBox(1,2) = MAXVAL(Mesh % Nodes % x)
        myBoundBox(2,1) = MINVAL(Mesh % Nodes % y)
        myBoundBox(2,2) = MAXVAL(Mesh % Nodes % y)
        myBoundBox(3,1) = MINVAL(Mesh % Nodes % z)
        myBoundBox(3,2) = MAXVAL(Mesh % Nodes % z)

        ! Wieso startet s bei 1????
        DO s = 1, nBodiesToComm
            DO i=0,totOFRanks-1
                ALLOCATE( OFp(i,s) % OFMesh )
                OFp(i,s) % globalRank = i + OFRanksStart
            END DO

            ALLOCATE( OF_EL_overlap(0:totLocalRanks-1,0:totOFRanks-1,nBodiesToComm) )
            ALLOCATE( ELboundBoxes(3,2,0:totLocalRanks-1,nBodiesToComm) )

            CALL findOverlappingBoxes(s)

            ! Starting communication
            !------------------------------------------------------------------------

            DO i = 0, totOFRanks - 1
                IF(.NOT.OFp(i,s) % boxOverlap) CYCLE
                ! Number of OpenFOAM cells
                CALL MPI_IRECV(OFp(i,s) % OFMesh % NumberOfNodes, 1, MPI_INTEGER, &
                                                OFp(i,s) % globalRank, 999, MPI_COMM_WORLD, OFp(i,s) % reqRecv, ierr)
            END DO

            DO i = 0, totOFRanks - 1
                IF(.NOT.OFp(i,s) % boxOverlap) CYCLE
                CALL MPI_TEST_SLEEP(OFp(i,s) % reqRecv, ierr)

                ALLOCATE( OFp(i,s) % OFMesh % Nodes, &
                                    OFp(i,s) % OFMesh % Variables, &
                                    OFp(i,s) % OFVar ) 
                ALLOCATE( OFp(i,s) % OFMesh % Nodes % x( OFp(i,s) % OFMesh % NumberOfNodes ), &
                                    OFp(i,s) % OFMesh % Nodes % y( OFp(i,s) % OFMesh % NumberOfNodes ), &
                                    OFp(i,s) % OFMesh % Nodes % z( OFp(i,s) % OFMesh % NumberOfNodes ), &
                                    OFp(i,s) % foundCells( OFp(i,s) % OFMesh % NumberOfNodes ) )

                OFp(i,s) % OFMesh % NumberOfBulkElements = 0
                OFp(i,s) % OFMesh % NumberOfBoundaryElements = 0
                OFp(i,s) % OFMesh % Projector => NULL()
                OFp(i,s) % foundCells = .FALSE.

                ! Cell x coordinates
                CALL MPI_IRECV(OFp(i,s) % OFMesh % Nodes % x, OFp(i,s) % OFMesh % NumberOfNodes, MPI_DOUBLE, &
                                             OFp(i,s) % globalRank, 997, MPI_COMM_WORLD, OFp(i,s) % reqRecv, ierr)
                CALL MPI_REQUEST_FREE(OFp(i,s) % reqRecv, ierr)
                CALL MPI_IRECV(OFp(i,s) % OFMesh % Nodes % y, OFp(i,s) % OFMesh % NumberOfNodes, MPI_DOUBLE, &
                                             OFp(i,s) % globalRank, 997, MPI_COMM_WORLD, OFp(i,s) % reqRecv, ierr)
                CALL MPI_REQUEST_FREE(OFp(i,s) % reqRecv, ierr)
                CALL MPI_IRECV(OFp(i,s) % OFMesh % Nodes % z, OFp(i,s) % OFMesh % NumberOfNodes, MPI_DOUBLE, &
                                             OFp(i,s) % globalRank, 997, MPI_COMM_WORLD, OFp(i,s) % reqRecv, ierr)
            END DO

            CALL Info('Elmer2OpenFOAMSolver','Projecting field to OpenFOAM cell centers',Level=10) 
            DO i = 0, totOFRanks - 1
                OFp(i,s) % nFoundCells = 0 ! keep this
                IF(.NOT.OFp(i,s) % boxOverlap) CYCLE
                ! wait for z coordinates
                CALL MPI_TEST_SLEEP(OFp(i,s) % reqRecv, ierr)

                IF ( CoordinateSystemDimension() == 2 ) THEN
                    IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
                             CurrentCoordinateSystem() == CylindricSymmetric ) THEN
                        OFp(i,s) % OFMesh % Nodes % x = SQRT(OFp(i,s) % OFMesh % Nodes % x**2 + OFp(i,s) % OFMesh % Nodes % z**2)
                    END IF
                    OFp(i,s) % OFMesh % Nodes % z = 0
                END IF

                CALL InterpolateMeshToMeshQ( OldMesh = Mesh, &
                        NewMesh                    = OFp(i,s) % OFMesh, &
                        UseQuadrantTree    = .TRUE., &
                        Projector                = OFp(i,s) % OFMesh % Projector, &
                        FoundNodes             = OFp(i,s) % foundCells, &
                        KeepUnfoundNodes = .FALSE.)

                OFp(i,s) % nFoundCells = COUNT(OFp(i,s) % foundCells)

                ! Number of cells found in each Elmer process
                CALL MPI_ISEND( OFp(i,s) % nFoundCells, 1, MPI_INTEGER, &
                                                OFp(i,s) % globalRank, 995, MPI_COMM_WORLD, OFp(i,s) % reqSend, ierr)
            END DO

            DO i = 0, totOFRanks - 1
                IF(.NOT.OFp(i,s) % boxOverlap) CYCLE
                ! wait for nFoundCells
                CALL MPI_TEST_SLEEP(OFp(i,s) % reqSend, ierr)

                IF ( OFp(i,s) % nFoundCells == 0 ) CYCLE
                ALLOCATE( OFp(i,s) % OFVar % Values( OFp(i,s) % nFoundCells ), &
                                    OFp(i,s) % foundCellsIndx( OFp(i,s) % nFoundCells ), &
                                    OFp(i,s) % OFVar % Perm( OFp(i,s) % nFoundCells ) )

                OFp(i,s) % OFVar % Perm = (/ (j, j = 1, OFp(i,s) % nFoundCells) /)
                OFp(i,s) % foundCellsIndx = PACK((/ (j, j = 0, OFp(i,s) % OFMesh % NumberOfNodes-1) /),OFp(i,s) % foundCells)

                ! Indexes for cells that were found on this piece of Elmer mesh
                CALL MPI_ISEND( OFp(i,s) % foundCellsIndx, OFp(i,s) % nFoundCells, MPI_INTEGER, &
                                                OFp(i,s) % globalRank, 994, MPI_COMM_WORLD, OFp(i,s) % reqSend, ierr)
            END DO

            DO i = 0, totOFRanks - 1
                IF ( OFp(i,s) % nFoundCells == 0 ) CYCLE
                CALL MPI_TEST_SLEEP(OFp(i,s) % reqSend, ierr)
            END DO
        END DO ! nBodiesToComm
    END IF ! .NOT. VISITED -> ran once


    ! run at com
    DO s = 1, nBodiesToComm
        ! Receive simulation status
        CALL Info('Elmer2OpenFOAM','Sending Body '//TRIM(I2S(s))//' Variables', Level=3 )

        CALL MPI_IRECV( OFstatus, 1, MPI_INTEGER, OFp(0,s) % globalRank, 799, MPI_COMM_WORLD, OFp(0,s) % reqRecv, ierr)
        CALL MPI_TEST_SLEEP(OFp(0,s) % reqRecv, ierr)

        IF (OFstatus.NE.1) THEN
            CALL Info('Elmer2OpenFOAM','Elmer has last iteration!', Level=3 )
            exitcond = ListGetCReal( CurrentModel % Simulation,'Exit Condition',Found)
                IF(.NOT.Found) CALL ListAddConstReal(CurrentModel % Simulation,'Exit Condition',1.0_dp)
        END IF

        ! Send fields
        DO j=1,nVars
            IF (.NOT. USE_TARGET_VARS_IN_BODY(j,s) .EQ. 0) THEN
                VarName = ListGetString( Params, 'Target Variable '//TRIM(I2S(j)), Found )
                Var => VariableGet( Mesh % Variables, VarName )

                IF(.NOT. ASSOCIATED( Var ) ) THEN
                    CALL Fatal('Elmer2OpenFOAMSolver','Variable '//TRIM(VarName)//' does not exist in Elmer mesh!')
                END IF

                DO i = 0, totOFRanks - 1
                    IF ( OFp(i,s) % nFoundCells == 0 ) CYCLE
                    OFp(i,s) % OFVar % Values = 0
                    ! value interpolation using fixed from weights?
                    CALL CRS_ApplyProjector( OFp(i,s) % OFMesh % Projector % Matrix, Var % Values, &
                                            Var % Perm, OFp(i,s) % OFVar % Values, OFp(i,s) % OFVar % Perm )

                    CALL MPI_ISEND( OFp(i,s) % OFVar % Values, OFp(i,s) % nFoundCells, MPI_DOUBLE, &
                                                    OFp(i,s) % globalRank, 1000, MPI_COMM_WORLD, OFp(i,s) % reqSend, ierr)
                END DO

                DO i = 0, totOFRanks - 1
                    IF ( OFp(i,s) % nFoundCells == 0 ) CYCLE
                    CALL MPI_TEST_SLEEP(OFp(i,s) % reqSend, ierr)
                END DO
            END IF
        END DO
    END DO ! nBodiesToComm

    VISITED = .TRUE.
    write(timeStr , '(F9.5)') MPI_WTIME() - commTime

    CALL Info('Elmer2OpenFOAM',' = '//TRIM(timeStr)//' s', Level=3 )

    CALL Info('Elmer2OpenFOAMSolver','All done', Level=4 )
    CALL Info('Elmer2OpenFOAMSolver','-----------------------------------------', Level=4 )
    
END SUBROUTINE Elmer2OpenFOAMSolver
