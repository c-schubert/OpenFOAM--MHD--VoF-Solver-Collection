if(updateElmerSolution || !runTime.run()) 
{
    // Receive fields form Elmer
    receiving.sendStatus(runTime.run());
    receiving.recvScalar(JH_recv);
}

