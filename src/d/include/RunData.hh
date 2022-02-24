//

/// \file RunData.hh
/// \brief Definition of the RunData class

#ifndef RunData_h
#define RunData_h 1

#include "G4Run.hh"
#include "globals.hh"

#include <array>
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int kDim = 36;

///  Run data class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers.
///
/// In order to reduce the number of data members a 2-dimensions array
/// is introduced for each quantity:
/// - fEdep[], fTrackLength[].
///
/// The data are collected step by step in SteppingAction, and
/// the accumulated values are filled in histograms and entuple
/// event by event in B4EventAction.

class RunData : public G4Run
{
public:
    RunData();
    virtual ~RunData();

    void SetTotDoseinfo(G4double accdosex, G4double acc_eq_dosex, G4double accphotonx);
    void SetTotParticleinfo(G4double elecnumx, G4double gammanumx, G4double protonnumx, G4double neutroncnumx, G4double alphanumx);
    void SetTotSipminfo(G4double e1, G4double e2, G4double c1, G4double c2);

    void SetEntPointG(G4double dx, G4double dy, G4double dz);
    void SetEntPointE(G4double dx, G4double dy, G4double dz);
    void SetExitPointG(G4double dx, G4double dy, G4double dz);
    void SetExitPointE(G4double dx, G4double dy, G4double dz);
    void SetEntAngle(G4double theta, G4double phi);
    void SetExitAngle(G4double theta, G4double phi);

    void detectorData1(G4double dx, G4double dy);
    void detectorData2(G4double dz);
    void FillPerEvent();
    void Reset();

    inline void sipm(G4double energy, G4int number)
    {
        G4double num = G4UniformRand();
        if ((number == 0 || number == 1) & num < 0.27)
        {
            fen[number] += energy;
            fcount[number]++;
        }
        num = G4UniformRand();
        if ((number == 2 || number == 3) & num < 0.50)
        {
            fen[number] += energy;
            fcount[number]++;
        }
    };

    //      G4cout<<"count"<<fcount[number]<<"  "<<tot_sipm_count[number]<<G4endl;

    inline void SetDepEnergy(G4double depen)
    {
        eDep += depen;
    }
    inline void SetDepEnergySmall(G4double depen)
    {
        eDep_s += depen;
    }

    inline void SetNumPhoton()
    {
        fcount_photon++;
    }

    inline void SetDose(G4double fdose)
    {
        dose += fdose;
    }
    inline void SetDoseSmall(G4double fdose)
    {
        dose_s += fdose;
    }
    inline void SetElectronenergy(G4double eenergy)
    {
        elecspec = eenergy;
        elecnum++;
    }
    inline void SetGammaenergy(G4double gerergy)
    {
        gammaspec = gerergy;
        gammanum++;
    }
    inline void SetProtonenergy(G4double perergy)
    {
        protonspec = perergy;
        protonnum++;
    }
    inline void SetNeutronenergy(G4double nerergy)
    {
        neutronspec = nerergy;
        neutronnum++;
    }
    inline void SetAlphaenergy(G4double aerergy)
    {
        alphaspec = aerergy;
        alphanum++;
    }
    inline void SetTotenergy(G4double toterergy)
    {
        totspec = toterergy;
    }
    inline G4double GetTotenergy()
    {
        return totspec;
    }
    inline void SetInitenergy(G4double ee)
    {
        initE = ee;
    }

    inline void SetGammainWallenergy(G4double nerergy)
    {
        GammainWallnum++;
    }
    inline void SetBetainWallenergy(G4double nerergy)
    {

        Betainwallnum++;
    }

    //////////////////////////////////////////////////////////////////////
    inline void SetEntElectronenergy(G4double eenergy)
    {
        entelecspec = eenergy;
        entelecnum++;
    }
    inline void SetEntGammaenergy(G4double gerergy)
    {
        entgammaspec = gerergy;
        entgammanum++;
    }

    inline void SetEntElectronenergySmall(G4double eenergy)
    {
        entelecnum_s++;
    }
    inline void SetEntGammaenergySmall(G4double gerergy)
    {
        entgammanum_s++;
    }
    inline void SetEntProtonenergy(G4double perergy)
    {
        entprotonspec = perergy;
        entprotonnum++;
    }
    inline void SetEntNeutronenergy(G4double nerergy)
    {
        entneutronspec = nerergy;
        entneutronnum++;
    }
    inline void SetEntAlphaenergy(G4double aerergy)
    {
        entalphaspec = aerergy;
        entalphanum++;
    }
    inline void SetEntTotenergy(G4double toterergy)
    {
        enttotspec = toterergy;
    }
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    inline void SetExitElectronenergy(G4double eenergy)
    {
        extelecspec = eenergy;
        extelecnum++;
    }
    inline void SetExitGammaenergy(G4double gerergy)
    {
        extgammaspec = gerergy;
        extgammanum++;
    }
    inline void SetExitProtonenergy(G4double perergy)
    {
        extprotonspec = perergy;
        extprotonnum++;
    }
    inline void SetExitNeutronenergy(G4double nerergy)
    {
        extneutronspec = nerergy;
        extneutronnum++;
    }
    inline void SetExitAlphaenergy(G4double aerergy)
    {
        extalphaspec = aerergy;
        extalphanum++;
    }
    inline void SetExitTotenergy(G4double toterergy)
    {
        exttotspec = toterergy;
    }
    //////////////////////////////////////////////////////
    inline void SetEqDose(G4double eqdose)
    {
        eq_dose += eqdose;
    }
    inline void SetEqDoseSmall(G4double eqdose)
    {
        eq_dose_s += eqdose;
    }

    inline G4int setCoincoef(G4int coef)
    {
        return coef;
    }

    inline G4int totNumEvent()
    {

        return (G4Run::GetNumberOfEventToBeProcessed());
    }

private:
    std::array<G4String, kDim> fVolumeNames;
    std::array<G4double, kDim> fEkin;
    std::array<G4double, kDim> fWavelength;

    G4double fen[4];
    G4int fcount[4];

    G4double fposx_exiG;
    G4double fposy_exiG;
    G4double fposz_exiG;
    G4double fposx_exiE;
    G4double fposy_exiE;
    G4double fposz_exiE;
    G4double fposx_entG;
    G4double fposy_entG;
    G4double fposz_entG;
    G4double fposx_entE;
    G4double fposy_entE;
    G4double fposz_entE;
    G4double phi_ent;
    G4double theta_ent;
    G4double phi_exi;
    G4double theta_exi;
    G4double eDep;
    G4double dose;
    G4double eDep_s;
    G4double dose_s;
    G4int fcount_photon;

    G4double eq_dose;
    G4double eq_dose_s;
    G4double elecspec;
    G4double protonspec;
    G4double gammaspec;
    G4double neutronspec;
    G4double alphaspec;
    G4double totspec; //
                      ///////////////////////////////////
    G4int elecnum;
    G4int protonnum;
    G4int gammanum;
    G4int neutronnum;
    G4int alphanum;
    ////////////////////////////////
    G4double entelecspec;
    G4double entprotonspec;
    G4double entgammaspec;
    G4double entneutronspec;
    G4double entalphaspec;
    G4double enttotspec; //
                         ///
                         ////////////////////////////////
    G4int entelecnum;
    G4int entelecnum_s;
    G4int entprotonnum;
    G4int entgammanum;
    G4int entgammanum_s;
    G4int entneutronnum;
    G4int entalphanum;
    /////////////////////////
    ////////////////////////////////
    G4double extelecspec;
    G4double extprotonspec;
    G4double extgammaspec;
    G4double extneutronspec;
    G4double extalphaspec;
    G4double exttotspec; //
                         ///
                         ////////////////////////////////
    G4int extelecnum;
    G4int extprotonnum;
    G4int extgammanum;
    G4int extneutronnum;
    G4int extalphanum;
    /////////////////////////
    G4int Betainwallnum;
    G4int GammainWallnum;
    G4double initE;
    //////////////////

    G4double accdose;
    G4double acc_eq_dose;
    G4int accphoton;
    G4int accelecnum;
    G4int accprotonnum;
    G4int accgammanum;
    G4int accneutronnum;
    G4int accalphanum;
    G4int tot_sipmc1;
    G4double tot_sipme1;
    G4int tot_sipmc2;
    G4double tot_sipme2;
    // std::vector<double*> tst;
};

// inline functions

// inline  void RunData::SetTotDoseinfo(G4double accdosex, G4double acc_eq_dosex, G4double accphotonx) {
//     accdose = accdosex;
//     acc_eq_dose = acc_eq_dosex;
//     accphoton = accphotonx;
// }
// inline  void RunData::SetTotParticleinfo(G4double elecnumx, G4double gammanumx
// , G4double protonnumx,G4double neutronnumx, G4double alphanumx) {
//     elecnum = elecnumx;
//     gammanum = gammanumx;
//     protonnum = protonnumx;
//     neutronnum = neutronnumx;
//     alphanum = alphanumx;
// }
// inline  void RunData::SetTotSipminfo(G4double e1, G4double e2, G4double c1, G4double c2) {
//     tot_sipmc1 = c1;
//     tot_sipmc2 = c2;
//     tot_sipme1 = e1;
//     tot_sipme2 = e2;
// }

inline void RunData::SetExitPointG(G4double dx, G4double dy, G4double dz)
{
    fposx_exiG = dx;
    fposy_exiG = dy;
    fposz_exiG = dz;
}
inline void RunData::SetExitPointE(G4double dx, G4double dy, G4double dz)
{
    fposx_exiE = dx;
    fposy_exiE = dy;
    fposz_exiE = dz;
}
inline void RunData::SetEntPointG(G4double dx, G4double dy, G4double dz)
{
    fposx_entG = dx;
    fposy_entG = dy;
    fposz_entG = dz;
}
inline void RunData::SetEntPointE(G4double dx, G4double dy, G4double dz)
{
    fposx_entE = dx;
    fposy_entE = dy;
    fposz_entE = dz;
}
inline void RunData::SetEntAngle(G4double theta, G4double phi)
{
    theta_ent = theta;
    phi_ent = phi;
}
inline void RunData::SetExitAngle(G4double theta, G4double phi)
{
    theta_exi = theta;
    phi_exi = phi;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif