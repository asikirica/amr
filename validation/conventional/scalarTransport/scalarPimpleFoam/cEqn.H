{
    volScalarField Dt ("Dt", (turbulence->nut()/ScT));

    for (label i = 0; i <= 4; i++)
    {
        fvScalarMatrix cEqn
        (
              fvm::ddt(c)
            + fvm::div(phi, c)
            - fvm::laplacian(Dt, c)
        );
        cEqn.relax();
        cEqn.solve();
    }
}
