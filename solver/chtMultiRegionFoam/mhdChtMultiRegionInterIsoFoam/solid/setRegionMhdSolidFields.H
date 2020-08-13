Elmer<fvMesh>& receiving  = solidReceiving[i];

volScalarField& JH_recv = solidJH_recv[i];
volScalarField& JH = solidJH[i];

JH.field() = JH_recv.field();
