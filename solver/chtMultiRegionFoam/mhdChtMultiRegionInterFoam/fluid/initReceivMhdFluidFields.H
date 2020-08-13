
forAll(fluidRegions, i)
{
    fluidReceiving[i].sendStatus(1);
    Info << "OF receiving fluid " << i << " JxB_recv from Elmer" << endl;
    fluidReceiving[i].recvVector(fluidJxB_recv[i]);
    Info << "OF receiving fluid " << i << " JH_recv from Elmer" << endl;
    fluidReceiving[i].recvScalar(fluidJH_recv[i]);

}
