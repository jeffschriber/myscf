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

extern "C"
PsiReturnType myscf(Options &options)
{
    int print = options.get_int("PRINT");
    int doTei = options.get_bool("DO_TEI");

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();

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

    // Form h = T + V by first cloning T and then adding V
    hMat->copy(tMat);
    hMat->add(vMat);


    int ntri = hMat->nrow();
    ntri = ntri*ntri;

    //define Vector of two electron integrals
    SharedVector TEI(new Vector("TEI array", ntri*(ntri+1.0)/2.0));

    if(doTei){
      //  psi::outfile->Printf("\n  Two-electron Integrals\n\n");

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
                    //psi::outfile->Printf("\t(%2d %2d | %2d %2d) = %20.15f\n",
                      //  p, q, r, s, buffer[intIter.index()]);
                    ++count;

                    int pq = INDEX2(p,q);
                    int rs = INDEX2(r,s);
                    int pqrs = INDEX2(pq,rs);
                    TEI->set(pqrs, buffer[intIter.index()]);
                }
            }
        }
        psi::outfile->Printf("\n\tThere are %d unique integrals\n\n", count);
    }
    //Diagonalize the AO-overlap matrix sMat
    int dim = sMat->nrow();

    boost::shared_ptr<Matrix> L(new Matrix("L",dim,dim));
    boost::shared_ptr<Matrix> sEvals(new Matrix("sEvals", dim, dim));
    boost::shared_ptr<Vector> lambda(new Vector("lambda", dim));

    sMat->diagonalize(L, lambda);
    sEvals->set_diagonal(lambda);
    sEvals->invert();
    sEvals->copy(sEvals->partial_cholesky_factorize(0.0,false)); //gets square root

    //get sMat^-1/2, called sMatInv
    boost::shared_ptr<Matrix> AB(new Matrix("AB", dim, dim));
    boost::shared_ptr<Matrix> sMatInv(new Matrix("sMatInv", dim, dim));

    AB->gemm(false, false, 1.0, L, sEvals, 0.0);
    sMatInv->gemm(false, true, 1.0, AB, L, 0.0);

    //Build the initial Fock Matrix
    boost::shared_ptr<Matrix> fInit(new Matrix("fInit", dim, dim));
    boost::shared_ptr<Matrix> fInitLeft(new Matrix("fInitLeft", dim, dim));

    fInitLeft->gemm(true, false, 1.0, sMatInv, hMat, 0);
    fInit->gemm(false, false, 1.0, fInitLeft, sMatInv, 0);

    //MOi is the initial Matrix resultant from Fock Matrix diagonalization, not
    //transformed to the original AO basis
    boost::shared_ptr<Matrix> MOi(new Matrix("MOi", dim, dim));
    boost::shared_ptr<Vector> initOrbEn(new Vector("initOrbEn", dim));

    fInit->diagonalize(MOi, initOrbEn);

    //Now, transform MOi to original AO basis
    boost::shared_ptr<Matrix> moCoefInit(new Matrix("moCoefInit", dim, dim));
    moCoefInit->gemm(false, false, 1.0, sMatInv, MOi, 0);

    //Build density matrix, summing over occupied MOs
    boost::shared_ptr<Matrix> dmInit(new Matrix("dmInit", dim, dim));
    int nocc = wfn->nalpha();

    dmInit->gemm('n','t', dim, dim, nocc, 1.0, moCoefInit, dim, moCoefInit, dim,0.0,dim);

    //Compute Initial SCF energy
    double Eelec = 0.0;
    double Einit;

    fInit->add(hMat);
    Eelec = dmInit->vector_dot(fInit);

    psi::outfile->Printf("\nInitial Electronic Energy: %8.12f a.u.\n", Eelec);

    Einit = Eelec + nucrep;
    psi::outfile->Printf("Initial HF-SCF Energy: %8.12f a.u.\n", Einit);

    //Initialize a bunch of matrices outside of loop

    boost::shared_ptr<Matrix> fMatLeft(new Matrix ("fMatLeft", dim, dim));
    boost::shared_ptr<Matrix> fMat(new Matrix ("Fock Matrix", dim, dim));
    boost::shared_ptr<Matrix> fMatOrtho(new Matrix ("fMatOrtho", dim, dim));
    boost::shared_ptr<Matrix> dMat(new Matrix ("Density Matrix", dim, dim));
    boost::shared_ptr<Matrix> dMatPrev(new Matrix ("dMatPrev", dim, dim));
    boost::shared_ptr<Matrix> cMatD(new Matrix ("cMatD", dim, dim));
    boost::shared_ptr<Matrix> cMat(new Matrix ("MO Coefficient Matrix", dim, dim));
    boost::shared_ptr<Vector> orben(new Vector("Orbital Energies", dim));

    dMat->copy(dmInit);

    int iter = 1;
    double dE = 1;
    double rmsD = 1;
    double Etot = 0.0;

    /***Start SCF procedure***/

    outfile->Printf("\nStarting SCF Procedure:\n");
    outfile->Printf("\nIter        E(elec)              E(tot)               Delta(E)              RMS(D)\n");

    while(fabs(dE) > 1e-12 and fabs(rmsD) > 7e-12){

        fMat->copy(hMat);

        for(int i=0; i<dim; i++ ){
            for(int j=0; j<dim; j++){
                for(int k=0; k<dim; k++){
                    for(int l=0; l<dim; l++){
                        int ij = INDEX2(i,j);
                        int kl = INDEX2(k,l);
                        int ijkl = INDEX2(ij,kl);
                        int ik = INDEX2(i,k);
                        int jl = INDEX2(j,l);
                        int ikjl = INDEX2(ik,jl);

                        fMat->add(i,j, dMat->get(k,l) * ( 2*TEI->get(ijkl) - TEI->get(ikjl) ) );

                    }
                }
            }
        }

        fMatOrtho->copy(fMat);

        //Build new density matrix
        fMatLeft->zero();
        cMatD->zero();
        cMat->zero();
        orben->zero();

        fMatLeft->gemm(true, false, 1.0, sMatInv, fMat, 0);
        fMat->gemm(false, false, 1.0, fMatLeft, sMatInv, 0);

        fMat->diagonalize(cMatD, orben);
        cMat->gemm(false, false, 1.0, sMatInv, cMatD, 0);

        dMatPrev->copy(dMat);
        dMat->zero();

        //form Density matrix from Ca
        dMat->gemm('n','t', dim, dim, nocc, 1.0, cMat, dim, cMat, dim,0.0,dim);

        //Compute new SCF energy
        Eelec = 0.0;
        double Eprev = Etot;

        hMat->add(fMatOrtho);
        Eelec = dMat->vector_dot(hMat);
        hMat->subtract(fMatOrtho);

        Etot = Eelec + nucrep;

        //Test for Convergence

        dE = Etot - Eprev;
        rmsD = 0;

        dMat->subtract(dMatPrev);
        rmsD = dMat->vector_dot(dMat);

        dMat->add(dMatPrev);
        rmsD = sqrt(rmsD);

        outfile->Printf("%02d %20.12f %20.12f %20.12f %20.12f\n", iter, Eelec, Etot, dE, rmsD);

        ++iter;
        if(iter == 200){
            outfile->Printf("SCF not Converged\n");
            break;
        }

    }
    outfile->Printf("\nSCF converged in %d cycles.\n", iter-1);

    outfile->Printf("\nSCF Electronic Energy: %20.12f a.u.", Eelec);
    outfile->Printf("\nSCF Nuclear Energy:    %20.12f a.u.", nucrep);
    outfile->Printf("\nSCF Total Energy:      %20.12f a.u.\n", Etot);

    //Calculate one electron properties

    outfile->Printf("\n =>One Electron properties<= \n");

    boost::shared_ptr<OEProp> prop(new OEProp());
    prop->set_title("SCF");
    prop->add("DIPOLE");
    prop->add("MULLIKEN_CHARGES");
    prop->compute();
    return Success;
}

}} // End Namespaces
