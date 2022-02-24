/// \file RunData.cc
/// \brief Implementation of the RunData class

#include "RunData.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EventManager.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunData::RunData()
    : G4Run()
{
    initE = 0;
    Betainwallnum = 0;
    GammainWallnum = 0;
    fposx_exiG = 0;
    fposy_exiG = 0;
    fposz_exiG = 0;
    fposx_exiE = 0;
    fposy_exiE = 0;
    fposz_exiE = 0;
    fposx_entG = 0;
    fposy_entG = 0;
    fposz_entG = 0;
    fposx_entE = 0;
    fposy_entE = 0;
    fposz_entE = 0;
    phi_ent = -9;
    theta_ent = -9;
    phi_exi = -9;
    theta_exi = -9;

    eDep = 0;
    dose = 0;
    eq_dose = 0;
    eDep_s = 0;
    dose_s = 0;
    eq_dose_s = 0;
    fcount_photon = 0;
    ////
    elecspec = 0;
    gammaspec = 0;
    protonspec = 0;
    neutronspec = 0;
    alphaspec = 0;
    totspec = 0; //
                 ////////
    entelecspec = 0;
    entgammaspec = 0;
    entelecnum_s = 0;
    entgammanum_s = 0;
    entprotonspec = 0;
    entneutronspec = 0;
    entalphaspec = 0;
    enttotspec = 0; //
    entelecnum = 0;
    entgammanum = 0;
    entneutronnum = 0;
    entalphanum = 0;
    entprotonnum = 0;
    /////////////////
    ////////
    extelecspec = 0;
    extgammaspec = 0;
    extprotonspec = 0;
    extneutronspec = 0;
    extalphaspec = 0;
    exttotspec = 0; //
    extelecnum = 0;
    extgammanum = 0;
    extneutronnum = 0;
    extalphanum = 0;
    extprotonnum = 0;
    /////////////////
    elecnum = 0;
    gammanum = 0;
    neutronnum = 0;
    alphanum = 0;
    protonnum = 0;
    for (G4int i_sipm = 0; i_sipm < 4; i_sipm++)
    {
        fen[i_sipm] = 0;
        fcount[i_sipm] = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunData::~RunData()
{
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunData::FillPerEvent()
{
    //   G4cout
    // <<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()
    // <<G4endl;
    // G4int count = 0;

    auto analysisManager = G4AnalysisManager::Instance();

    //  G4GeneralParticleSource* fParticleGun =  new G4GeneralParticleSource();
    //     G4double initE=fParticleGun->GetParticleEnergy()/keV;

    //////////////////G4double initE=totspec;
    // G4cout<<initE<<G4endl;
    if (fcount[0] > 2)
    {
        analysisManager->FillH1(0, initE, 1);

        analysisManager->FillH1(1, initE, fen[0]);
    }
    if (fcount[1] > 2)
    {
        analysisManager->FillH1(2, initE, 1);

        analysisManager->FillH1(3, initE, fen[1]);
    }
    ////////
    if (fcount[2] > 2)
    {
        analysisManager->FillH1(4, initE, 1);

        analysisManager->FillH1(5, initE, fen[2]);
    }
    if (fcount[3] > 2)
    {
        analysisManager->FillH1(6, initE, 1);

        analysisManager->FillH1(7, initE, fen[3]);
    }
    ////////////////////////////////////////// coincidencebig
    if (fcount[0] > 2 && fcount[1] > 2)
    {
        analysisManager->FillH1(8, initE, 1);

        analysisManager->FillH1(9, initE, fen[0] + fen[1]);

        ///

        ///
        analysisManager->FillH1(10, initE, eDep);

        analysisManager->FillH1(11, initE, dose);
    }
    //////////////////////////////////////////////////////// coincidencesmall
    if (fcount[2] > 2 && fcount[3] > 2)
    {
        analysisManager->FillH1(12, initE, 1);

        analysisManager->FillH1(13, initE, fen[2] + fen[3]);

        ///
        ///
        analysisManager->FillH1(14, initE, eDep_s);

        analysisManager->FillH1(15, initE, dose_s);
    }
    //////
    if (dose > 0)
    {
        analysisManager->FillH1(16, initE, eDep);

        analysisManager->FillH1(17, initE, dose);

        analysisManager->FillH1(18, initE, eq_dose);
    }
    //////
    if (dose_s > 0)
    {
        analysisManager->FillH1(19, initE, eDep_s);

        analysisManager->FillH1(20, initE, dose_s);

        analysisManager->FillH1(21, initE, eq_dose_s);
    }
    //////////////////////////////////////////gamma
    if (entgammanum != 0)
    {
        analysisManager->FillH1(22, initE, 1);
    }
    if (entelecnum != 0)
    {
        analysisManager->FillH1(23, initE, 1);
    }
    if (entgammanum_s != 0)
    {
        analysisManager->FillH1(24, initE, 1);
    }
    if (entelecnum_s != 0)
    {
        analysisManager->FillH1(25, initE, 1);
    }
    /////////////////////////////////////
    if (GammainWallnum != 0)
    {
        analysisManager->FillH1(26, initE, 1);
    }
    if (Betainwallnum != 0)
    {
        analysisManager->FillH1(27, initE, 1);
    }

    analysisManager->FillH1(28, initE, 1);

    // }
    // /////
    // if (entgammanum!=0)
    // {
    // analysisManager->FillH1(count,initE, 1);
    //     count++;
    // /////traveling!!!!!
    // // analysisManager->FillH1(count,initE, sqrt(pow(fposx_exiG-fposx_entG,2)+pow(fposy_exiG-fposy_entG,2)+pow(fposz_exiG-fposz_entG,2)));
    // //     count++;
    // //G4cout<<sqrt(pow(fposx_exiG-fposx_entG,2)+pow(fposy_exiG-fposy_entG,2)+pow(fposz_exiG-fposz_entG,2))<<G4endl;
    // }

    // /////////////////////////////////////////////// electron
    // if (dose>0)
    // {
    // if (entelecnum!=0)
    // {
    // analysisManager->FillH1(count,initE, 1);
    //     count++;
    // }
    // }
    // if (entelecnum!=0)
    // {
    // analysisManager->FillH1(count,initE, 1);
    //     count++;
    // /////traveling!!!!!
    // // analysisManager->FillH1(count,initE, sqrt(pow(fposx_exiE-fposx_entE,2)+pow(fposy_exiE-fposy_entE,2)+pow(fposz_exiE-fposz_entE,2)));
    // //     count++;
    // }

    // analysisManager->FillH1(count,initE, 1);
    //     count++;

    // accumulate statistic
    // in the order od the histograms, ntuple columns declarations
    //                     if (fcount[0]>0)  {
    //    analysisManager->FillNtupleDColumn(counter++,fcount[0]);
    //     analysisManager->FillNtupleDColumn(counter++,fen[0]);
    //                     }
    //                     else{
    //    analysisManager->FillNtupleDColumn(counter++,fcount[0]);
    //     analysisManager->FillNtupleDColumn(counter++,0);
    //                     }
    ////////////////////////////////////////////////////////////////////
    //                      if (fcount[1]>0)  {
    //    analysisManager->FillNtupleDColumn(counter++,fcount[1]);
    //     analysisManager->FillNtupleDColumn(counter++,fen[1]);
    //                     }
    //                     else{
    //    analysisManager->FillNtupleDColumn(counter++,0/*fcount[0]*/);
    //     analysisManager->FillNtupleDColumn(counter++,0);
    //                     }
    //                 if (fcount[0]>2&&fcount[1]>2)  {
    //     for (G4int i_sipm = 0; i_sipm < 2; i_sipm++) {
    //    analysisManager->FillNtupleDColumn(counter++,fcount[i_sipm]);
    //     analysisManager->FillNtupleDColumn(counter++,fen[i_sipm]);

    //     }
    //                 }
    //                  else {
    //           for (G4int i_sipm = 0; i_sipm < 2; i_sipm++) {
    //         analysisManager->FillNtupleDColumn(counter++,0);
    //             analysisManager->FillNtupleDColumn(counter++,0);
    //              }
    //                  }

    //     analysisManager->FillNtupleDColumn(counter++,fposx_ent);
    //     analysisManager->FillNtupleDColumn(counter++,fposy_ent);
    //     analysisManager->FillNtupleDColumn(counter++,fposz_ent);
    //     ////
    //     analysisManager->FillNtupleDColumn(counter++,fposx_exi);
    //     analysisManager->FillNtupleDColumn(counter++,fposy_exi);
    //     analysisManager->FillNtupleDColumn(counter++,fposz_exi);
    //     //
    //     analysisManager->FillNtupleDColumn(counter++,fposx_sci);
    //     analysisManager->FillNtupleDColumn(counter++,fposy_sci);
    //     analysisManager->FillNtupleDColumn(counter++,fposz_sci);
    //     ///
    //         analysisManager->FillNtupleDColumn(counter++,theta_ent);
    //     analysisManager->FillNtupleDColumn(counter++,phi_ent);
    //         analysisManager->FillNtupleDColumn(counter++,theta_exi);
    //         analysisManager->FillNtupleDColumn(counter++,phi_exi);
    //         ////
    //     analysisManager->FillNtupleDColumn(counter++,eDep);
    //         analysisManager->FillNtupleDColumn(counter++,dose);
    //               analysisManager->FillNtupleDColumn(counter++,eq_dose);
    // //
    //     analysisManager->FillNtupleDColumn(counter++,elecspec);
    //     analysisManager->FillNtupleDColumn(counter++,gammaspec);
    //             analysisManager->FillNtupleDColumn(counter++,alphaspec);
    //                     analysisManager->FillNtupleDColumn(counter++,protonspec);
    //                 analysisManager->FillNtupleDColumn(counter++,neutronspec);
    //                 analysisManager->FillNtupleDColumn(counter++,totspec);
    // //
    //  analysisManager->FillNtupleDColumn(counter++,entelecspec);
    //     analysisManager->FillNtupleDColumn(counter++,entgammaspec);
    //             analysisManager->FillNtupleDColumn(counter++,entalphaspec);
    //                     analysisManager->FillNtupleDColumn(counter++,entprotonspec);
    //                 analysisManager->FillNtupleDColumn(counter++,entneutronspec);
    //                 analysisManager->FillNtupleDColumn(counter++,enttotspec);
    // //
    // //
    //  analysisManager->FillNtupleDColumn(counter++,extelecspec);
    //     analysisManager->FillNtupleDColumn(counter++,extgammaspec);
    //             analysisManager->FillNtupleDColumn(counter++,extalphaspec);
    //                     analysisManager->FillNtupleDColumn(counter++,extprotonspec);
    //                 analysisManager->FillNtupleDColumn(counter++,extneutronspec);
    //                 analysisManager->FillNtupleDColumn(counter++,exttotspec);
    //  //

    //                 //
    //                   analysisManager->FillNtupleDColumn(counter++,elecnum);
    //                  analysisManager->FillNtupleDColumn(counter++,gammanum);
    //                analysisManager->FillNtupleDColumn(counter++,neutronnum);
    //                 analysisManager->FillNtupleDColumn(counter++,protonnum);
    //                 analysisManager->FillNtupleDColumn(counter++,alphanum);
    //  //

    //                 //
    //                   analysisManager->FillNtupleDColumn(counter++,entelecnum);
    //                  analysisManager->FillNtupleDColumn(counter++,entgammanum);
    //                analysisManager->FillNtupleDColumn(counter++,entneutronnum);
    //                 analysisManager->FillNtupleDColumn(counter++,entprotonnum);
    //                 analysisManager->FillNtupleDColumn(counter++,entalphanum);
    //  //

    //                 //
    //                   analysisManager->FillNtupleDColumn(counter++,extelecnum);
    //                  analysisManager->FillNtupleDColumn(counter++,extgammanum);
    //                analysisManager->FillNtupleDColumn(counter++,extneutronnum);
    //                 analysisManager->FillNtupleDColumn(counter++,extprotonnum);
    //                 analysisManager->FillNtupleDColumn(counter++,extalphanum);
    //                 analysisManager->FillNtupleDColumn(counter++,fcount_photon);

    // //             if (fcount[0]>0&&fcount[1]>0)  {
    // //  analysisManager->FillNtupleDColumn(counter++,1/*fcount[0]*//*+fcount[1]*/);////////////////////////////////////////////////////
    // //   analysisManager->FillNtupleDColumn(counter++,fen[0]/*+fen[1])/1000*/);////////////////////
    // //             }
    // //              else {  analysisManager->FillNtupleDColumn(counter++,0);
    // //             analysisManager->FillNtupleDColumn(counter++,0);
    // //              }

    //     analysisManager->AddNtupleRow();

    Reset();

    // //G4cout<<num_Cphotons_ck<<"  "<<num_Cphotons<<G4endl;

    //   if ((G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()+1)==totNumEvent())

    //   {
    //      G4cout<<"- Total dose(Gy)    Dose equivalent rate (Sv)  Dose equivalent rate (Sv/h)-----------"<<G4endl;
    //       G4cout<<accdose<<" "<<acc_eq_dose<<" "<<acc_eq_dose*3600<<G4endl;
    //        G4cout<<"---------------------------------total photons created-------------------------------"<<G4endl;
    //      G4cout<<accphoton<<G4endl;
    //        G4cout<<"---------------------------------entered particles-------------------------------"<<G4endl;
    //     //  G4cout<<tot_fcoince<<G4endl;
    //         G4cout<<"- No.of incident e gamma alpha neutron proton------------------------------"<<G4endl;
    //       G4cout<<accelecnum<<" "<<accgammanum<<" "<<accalphanum<<" "<<accneutronnum<<" "<<accprotonnum<<G4endl;
    //         G4cout<<"---------------------------------total energy in sipm---------------------"<<G4endl;

    //         G4cout<<"In sipm1: "<<tot_sipme1/1e3<<" keV"<<G4endl;
    //         G4cout<<"In sipm2: "<<tot_sipme2/1e3<<" keV"<<G4endl;
    //        G4cout<<"---------------------------------cpm in sipm-------------------------------"<<G4endl;
    //         G4cout<<"In sipm1: "<<tot_sipmc1<<G4endl;
    //         G4cout<<"In sipm2: "<<tot_sipmc2<<G4endl;

    //   }
    // tot_fcoince ++;
    //  if (G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()==0){
    //       tot_dose=0;
    //     tot_photon=0;
    //     tot_eq_dose=0;
    //          accelecnum=0;
    //      accgammanum=0;
    //      accneutronnum=0;
    //      accalphanum=0;
    //           accprotonnum=0;
    //                tot_sipmcount1=0;
    //      tot_sipmenergy1=0;
    //        tot_sipmcount2=0;
    //     tot_sipmenergy2=0;
    //  }

    // G4cout<<tot_dose<<" "<<tot_photon<<" "<<tot_eq_dose<<""<<G4endl;
    //     //tot_dose+=dose;
    //   //  tot_eq_dose+=eq_dose;
    // //tot_photon+=fcount_photon;
    //   //  accelecnum+=elecnum;
    //    //  accprotonnum+=protonnum;
    //     //accgammanum+=gammanum;
    //     //accneutronnum+=neutronnum;
    //     //accalphanum+=alphanum;
    //  analysisManager->AddNtupleRow();

    //  if ((G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()+1)==
    //         G4Run::GetNumberOfEventToBeProcessed())
    //  {
    //      G4cout<<"- Total dose(Gy)    Dose equivalent rate (Sv)  Dose equivalent rate (Sv/h)-----------"<<G4endl;
    //       G4cout<<tot_dose<<" "<<tot_eq_dose<<" "<<tot_eq_dose*3600<<G4endl;
    //        G4cout<<"---------------------------------total photons created-------------------------------"<<G4endl;
    //      G4cout<<tot_photon<<G4endl;
    //        G4cout<<"---------------------------------coincidence coefficient in the events-------------------------------"<<G4endl;
    //      G4cout<<tot_fcoince<<G4endl;
    //         G4cout<<"- No.of incident e gamma alpha neutron proton------------------------------"<<G4endl;
    //       G4cout<<accelecnum<<" "<<accgammanum<<" "<<accalphanum<<" "<<accneutronnum<<" "<<accprotonnum<<G4endl;
    //         G4cout<<"---------------------------------total energy in sipm---------------------"<<G4endl;

    //         G4cout<<"In sipm1: "<<tot_sipmenergy1/1e3<<" keV"<<G4endl;
    //         G4cout<<"In sipm2: "<<tot_sipmenergy2/1e3<<" keV"<<G4endl;
    //        G4cout<<"---------------------------------cpm in sipm-------------------------------"<<G4endl;
    //         G4cout<<"In sipm1: "<<tot_sipmcount1<<" keV"<<G4endl;
    //         G4cout<<"In sipm2: "<<tot_sipmcount2<<" keV"<<G4endl;

    //  G4cout<<"- Total dose(Gy)   Dose equivalent rate (Sv/h)------------------------------"<<G4endl;
    // G4cout<<tot_dose<<" "<<tot_eq_dose*3600<<G4endl;
    //   G4cout<<"- No.of incident e gamma alpha neutron proton------------------------------"<<G4endl;
    //      G4cout<<elecnum<<" "<<gammanum<<" "<<alphanum<<" "<<neutronnum<<" "<<protonnum<<G4endl;
    //        G4cout<<"---------------------------------total energy in sipm---------------------"<<G4endl;
    //            for (G4int i_sipm = 0; i_sipm <(sizeof(tot_sipm_energy)/sizeof(*tot_sipm_energy)); i_sipm++) {
    //         G4cout<<"In sipm"<<i_sipm+1<<": "<<tot_sipm_energy[i_sipm]/1e3<<" keV"<<G4endl;}

    //        G4cout<<"---------------------------------cpm in sipm-------------------------------"<<G4endl;
    //            for (G4int i_sipm = 0; i_sipm <(sizeof(tot_sipm_count)/sizeof(*tot_sipm_count)); i_sipm++) {
    //         G4cout<<"In sipm"<<i_sipm+1<<": "<<tot_sipm_count[i_sipm]*60<<" cpm"<<G4endl;}
    //                G4cout<<"---------------------------------cpm in sipm-------------------------------"<<G4endl;
    //    G4cout<<"---------------------------------total photons created-------------------------------"<<G4endl;
    //  G4cout<<tot_photon<<G4endl;
    //    G4cout<<"---------------------------------coincidence coefficient in the events-------------------------------"<<G4endl;
    //  G4cout<<tot_fcoince<<G4endl;
    //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunData::Reset()
{
    for (G4int i_sipm = 0; i_sipm < 4; i_sipm++)
    {
        fen[i_sipm] = 0;
        fcount[i_sipm] = 0;
    }
    fposx_exiG = 0;
    fposy_exiG = 0;
    fposz_exiG = 0;
    fposx_exiE = 0;
    fposy_exiE = 0;
    fposz_exiE = 0;
    fposx_entG = 0;
    fposy_entG = 0;
    fposz_entG = 0;
    fposx_entE = 0;
    fposy_entE = 0;
    fposz_entE = 0;
    phi_ent = -9;
    theta_ent = -9;
    phi_exi = -9;
    theta_exi = -9;

    eDep = 0;
    dose = 0;
    eq_dose = 0;
    eDep_s = 0;
    dose_s = 0;
    eq_dose_s = 0;
    fcount_photon = 0;
    ////
    elecspec = 0;
    gammaspec = 0;
    protonspec = 0;
    neutronspec = 0;
    alphaspec = 0;
    totspec = 0; //
                 ////////
    entelecspec = 0;
    entgammaspec = 0;
    entprotonspec = 0;
    entneutronspec = 0;
    entalphaspec = 0;
    enttotspec = 0; //
    entelecnum = 0;
    entgammanum = 0;
    entelecnum_s = 0;
    entgammanum_s = 0;
    entneutronnum = 0;
    entalphanum = 0;
    entprotonnum = 0;
    /////////////////
    ////////
    extelecspec = 0;
    extgammaspec = 0;
    extprotonspec = 0;
    extneutronspec = 0;
    extalphaspec = 0;
    exttotspec = 0; //
    extelecnum = 0;
    extgammanum = 0;
    extneutronnum = 0;
    extalphanum = 0;
    extprotonnum = 0;
    /////////////////
    initE = 0;
    Betainwallnum = 0;
    GammainWallnum = 0;
    elecnum = 0;
    gammanum = 0;
    neutronnum = 0;
    alphanum = 0;
    protonnum = 0;
    for (G4int i_sipm = 0; i_sipm < 4; i_sipm++)
    {
        fen[i_sipm] = 0;
        fcount[i_sipm] = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
