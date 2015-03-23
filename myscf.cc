/*
 *@BEGIN LICENSE
 *
 * myscf by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>

INIT_PLUGIN

typedef boost::shared_ptr<psi::Matrix> SharedMatrix;

namespace psi{ namespace myscf {

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "MYSCF"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- Whether to compute two-electron integrals -*/
        options.add_bool("DO_TEI", true);
    }

    return true;
}

#define INDEX(i,j) (i>j) ? (ioff[i]+j):(ioff[j]+i)

extern "C"
PsiReturnType myscf(Options &options)
{
    int print = options.get_int("PRINT");
    int doTei = options.get_bool("DO_TEI");

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Form basis object:
    // Create a basis set parser object.
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");

    // The integral factory oversees the creation of integral objects
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory
            (aoBasis, aoBasis, aoBasis, aoBasis));
    
    // N.B. This should be called after the basis has been built, because the geometry has not been
    // fully initialized until this time.
    molecule->print();
    int nbf[] = { aoBasis->nbf() };
    double nucrep = molecule->nuclear_repulsion_energy();
    psi::outfile->Printf("\n    Nuclear repulsion energy: %16.8f\n\n", nucrep);

    // The matrix factory can create matrices of the correct dimensions...
    boost::shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(1, nbf, nbf);

    // Form the one-electron integral objects from the integral factory
    boost::shared_ptr<OneBodyAOInt> sOBI(integral->ao_overlap());
    boost::shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
    boost::shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(factory->create_matrix("Overlap"));
    SharedMatrix tMat(factory->create_matrix("Kinetic"));
    SharedMatrix vMat(factory->create_matrix("Potential"));
    SharedMatrix hMat(factory->create_matrix("One Electron Ints"));
    // Compute the one electron integrals, telling each object where to store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);

    sMat->print();
    tMat->print();
    vMat->print();

    // Form h = T + V by first cloning T and then adding V
    hMat->copy(tMat);
    hMat->add(vMat);
    hMat->print();

    int ntri = hMat->nrow();
    ntri = ntri*ntri;

    Vector ioff("ioff",ntri*(ntri+1.0)/2.0 );

    ioff[0]=0;
    for(int i=1; i<ntri*(ntri+1.0)/2.0; i++){
        ioff[i]=ioff[i-1]+i;
    }

    SharedVector TEI(new Vector("TEI array", ntri*(ntri+1.0)/2.0));

    if(doTei){
        psi::outfile->Printf("\n  Two-electron Integrals\n\n");

        // Now, the two-electron integrals
        boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
        // The buffer will hold the integrals for each shell, as they're computed
        const double *buffer = eri->buffer();
        // The iterator conveniently lets us iterate over functions within shells
        AOShellCombinationsIterator shellIter = integral->shells_iterator();
        int count=0;
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
            // Compute quartet
            if (eri->compute_shell(shellIter)) {
                // From the quartet get all the integrals
                AOIntegralsIterator intIter = shellIter.integrals_iterator();
                for (intIter.first(); intIter.is_done() == false; intIter.next()) {
                    int p = intIter.i();
                    int q = intIter.j();
                    int r = intIter.k();
                    int s = intIter.l();
                    psi::outfile->Printf("\t(%2d %2d | %2d %2d) = %20.15f\n",
                        p, q, r, s, buffer[intIter.index()]);
                    ++count;

                    int pq = INDEX(p,q);
                    int rs = INDEX(r,s);
                    int pqrs = INDEX(pq,rs);
                    TEI->set(pqrs, buffer[intIter.index()]);
                }
            }
        }
        psi::outfile->Printf("\n\tThere are %d unique integrals\n\n", count);
    }
    //Diagonalize the AO-overlap matrix sMat
    int dim = sMat->nrow();

    Matrix L("L", dim, dim);
    Matrix sEvals("sEvals", dim, dim);
    Vector lambda("lambda", dim);

    sMat->diagonalize(L, lambda);

    //Map eigenvalue vector onto a diagonal matrix
    for(int i=0; i<dim;i++){
        for(int j=0; j<dim; j++){
            if(i==j){
                sEvals.set(i, j, sqrt(lambda.get(j)) );
            }
            else{
                sEvals.set(i,j,0);
            }
        }
    }


    sEvals.invert();
    sEvals.print();

    //get sMat^-1/2, called sMatInv

    Matrix AB("AB", dim, dim);
    Matrix sMatInv("sMatInv", dim, dim);

    AB.gemm(false, false, 1.0, L, sEvals, 0);
    sMatInv.gemm(false, true, 1.0, AB, L, 0);
    sMatInv.print();


    //Build the initial Fock Matrix

    Matrix fInit("fInit", dim, dim);
    Matrix fInitLeft("fInitLeft", dim, dim);

    fInitLeft.gemm(true, false, 1.0, sMatInv, hMat, 0);
    fInit.gemm(false, false, 1.0, fInitLeft, sMatInv, 0);


    //MOi is the initial Matrix resultant from Fock Matrix diagonalization, not
    //transformed to the original AO basis

    Matrix MOi("MOi", dim, dim);
    Vector initOrbEn("initOrbEn", dim);

    fInit.diagonalize(MOi, initOrbEn);

    fInit.print();
    //Now, transform MOi to original AO basis

    Matrix moCoefInit("moCoefInit", dim, dim);

    moCoefInit.gemm(false, false, 1.0, sMatInv, MOi, 0);

    //Build density matrix using occupied MOs
    Matrix dmInit("dmInit", dim, dim);

    int natom = molecule->natom();
    int Ztot = 0;

    for(int i=0; i<natom; i++){
        Ztot += molecule->fZ(i);
    }

    int nocc = Ztot/2;
    psi::outfile->Printf("nocc: %d", nocc );

    for(int v=0; v<dim; v++){
        for(int u=0; u<dim; u++){
            for(int m=0; m<nocc;m++){
                    dmInit.add(u,v, moCoefInit(u,m) * moCoefInit(v,m) );
            }
        }
    }

    //Compute Initial SCF energy

    double Eelec = 0.0;
    double Einit;


    for(int u=0; u<dim; u++){
        for(int v=0; v<dim; v++){
            Eelec += dmInit.get(u,v)*(2.0*hMat->get(u,v) );
        }
    }

    psi::outfile->Printf("\nInitial Electronic Energy: %8.12f a.u.\n", Eelec);

    Einit = Eelec + nucrep;
    psi::outfile->Printf("Initial HF-SCF Energy: %8.12f a.u.\n", Einit);

    //Initialize a bunch of matrices outside of loop

    Matrix fMatLeft("fMatLeft", dim, dim);
    Matrix fMat("Fock Matrix", dim, dim);
    Matrix fMatOrtho("fMatOrtho", dim, dim);
    Matrix dMat("Density Matrix", dim, dim);
    Matrix dMatPrev("dMatPrev", dim, dim);
    Matrix cMatD("cMatD", dim, dim);
    Matrix cMat("MO Coefficient Matrix", dim, dim);
    Vector orben("Orbital Energies", dim);

    dMat.copy(dmInit);

    int iter = 1;
    double dE = 1;
    double rmsD = 1;
    double Etot;

    //Start SCF procedure

    outfile->Printf("\nStarting SCF Procedure:\n");
    outfile->Printf("\nIter        E(elec)              E(tot)               Delta(E)              RMS(D)\n");

        while(fabs(dE) > 1e-12 && fabs(rmsD) > 7e-12){

        fMat.copy(hMat);

        for(int i=0; i<dim; i++ ){
            for(int j=0; j<dim; j++){
                for(int k=0; k<dim; k++){
                    for(int l=0; l<dim; l++){
                        int ij = INDEX(i,j);
                        int kl = INDEX(k,l);
                        int ijkl = INDEX(ij,kl);
                        int ik = INDEX(i,k);
                        int jl = INDEX(j,l);
                        int ikjl = INDEX(ik,jl);

                        fMat.add(i,j, dMat.get(k,l) * ( 2*TEI->get(ijkl) - TEI->get(ikjl) ) );

                    }
                }
            }
        }

        fMatOrtho.copy(fMat);

        //Build new density matrix
        fMatLeft.zero();
        cMatD.zero();
        cMat.zero();
        orben.zero();

        fMatLeft.gemm(true, false, 1.0, sMatInv, fMat, 0);
        fMat.gemm(false, false, 1.0, fMatLeft, sMatInv, 0);

        fMat.diagonalize(cMatD, orben);
        cMat.gemm(false, false, 1.0, sMatInv, cMatD, 0);

        dMatPrev.copy(dMat);
        dMat.zero();
        for(int u=0; u<dim; u++){
            for(int v=0; v<dim; v++){
                for(int m=0; m<nocc; m++){
                    dMat.add(u,v, cMat.get(u,m)*cMat.get(v,m));
                    }
                }
            }


        //Compute new SCF energy
        if(iter ==1){
            Etot = Einit;
        }
        double Eprev = Etot;
        Etot = 0;
        Eelec = 0;

        for(int u=0; u<dim; u++){
            for(int v=0; v<dim; v++){
                Eelec += dMat.get(u,v) * ( hMat->get(u,v) + fMatOrtho.get(u,v) );
            }
        }

        Etot = Eelec + nucrep;

        //Test for Convergence

        dE = Etot - Eprev;
        rmsD = 0;
        for(int i=0; i<dim; i++){
            for(int j=0; j<dim; j++){
                rmsD += (dMat.get(i,j)-dMatPrev.get(i,j)) * (dMat.get(i,j)-dMatPrev.get(i,j));
            }
        }

        rmsD = sqrt(rmsD);


        outfile->Printf("%02d %20.12f %20.12f %20.12f %20.12f\n", iter, Eelec, Etot, dE, rmsD);

        ++iter;
        if(iter == 200){
            outfile->Printf("SCF not Converged\n");
            break;
        }

    }

    outfile->Printf("\nSCF Electronic Energy: %20.12f a.u.", Eelec);
    outfile->Printf("\nSCF Nuclear Energy:    %20.12f a.u.", nucrep);
    outfile->Printf("\nSCF Total Energy:      %20.12f a.u.\n", Etot);

    //Calculate on electron properties moment from components

    outfile->Printf("\n =>One Electron properties<= \n");

    boost::shared_ptr<OEProp> prop(new OEProp());
    prop->set_title("SCF");
    prop->add("DIPOLE");
    prop->add("MULLIKEN_CHARGES");
    prop->compute();
    return Success;
}

}} // End Namespaces
