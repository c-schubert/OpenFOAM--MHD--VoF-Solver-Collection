forAll(solidRegions, i)
{
    solidReceiving[i].sendStatus(1);
    Info << "OF receiving solid " << i << " JH_recv from Elmer" << endl;
    solidReceiving[i].recvScalar(solidJH_recv[i]);
}
