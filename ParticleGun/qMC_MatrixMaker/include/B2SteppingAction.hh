//
//  Header.h
//  
//
//  Created by Shiv Akshar Yadavalli on 15/02/16.
//
//

#ifndef _Header_h
#define _Header_h

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingVerbose.hh"


///////////////////////////
class B2SteppingAction : public G4UserSteppingAction
///////////////////////////
{
    
    //--------
public:
    //--------
    
    // Constructor and destructor
    B2SteppingAction(){}
    virtual ~ B2SteppingAction(){}
    
    // Member functions
    virtual void UserSteppingAction(const G4Step* fStep){
        G4Track * fTrack = fStep->GetTrack();
    G4int prec = G4cout.precision(16);
        
        G4cout<<G4endl;

        G4cout.precision(prec);
        
    
        
    }
    
    //-----------
protected:
    //-----------
    
    // Member data
    G4SteppingManager* fpSteppingManager;
    
};







#endif
