/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "myLinearSpring_polyParabolic.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(myLinearSpring_polyParabolic, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        myLinearSpring_polyParabolic,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::myLinearSpring_polyParabolic::myLinearSpring_polyParabolic
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict),
    anchor_(),
    refAttachmentPt_(),
    stiffness_(),
    damping_(),
    restLength_()
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::myLinearSpring_polyParabolic::~myLinearSpring_polyParabolic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::myLinearSpring_polyParabolic::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.transform(refAttachmentPt_);

    vector r = restraintPosition - anchor_;

    scalar magR = mag(r);
    r /= (magR + vSmall);
    scalar length_m = 1.045*restLength_;

    vector v = motion.velocity(restraintPosition);
    
    if ( magR > restLength_ && magR < length_m )
    {
	    scalar exp_z = Foam::pow((magR - restLength_) / (length_m - restLength_), 3);
	    restraintForce = -(10 + stiffness_ * exp_z)*(magR - restLength_)*r - damping_*(r & v)*r;
    }
	
    else

    {
	    restraintForce = (magR - magR)*r;
    }

    restraintMoment = Zero;

    if (motion.report())
    {
        Info<< " anchor point " << anchor_
            << " attachmentPt " << restraintPosition
            << " attachmentPt - anchor " << r*magR
            << " spring length " << magR
            << " force " << restraintForce
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::myLinearSpring_polyParabolic::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.lookup("anchor") >> anchor_;
    sDoFRBMRCoeffs_.lookup("refAttachmentPt") >> refAttachmentPt_;
    sDoFRBMRCoeffs_.lookup("stiffness") >> stiffness_;
    sDoFRBMRCoeffs_.lookup("damping") >> damping_;
    sDoFRBMRCoeffs_.lookup("restLength") >> restLength_;

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::myLinearSpring_polyParabolic::write
(
    Ostream& os
) const
{
    writeEntry(os, "anchor", anchor_);

    writeEntry(os, "refAttachmentPt", refAttachmentPt_);

    writeEntry(os, "stiffness", stiffness_);

    writeEntry(os, "damping", damping_);

    writeEntry(os, "restLength", restLength_);
}

// ************************************************************************* //
