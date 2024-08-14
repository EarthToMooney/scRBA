**** OUTPUT FILES ***********************************************
*   FVA.txt -        GAMS output text file that         *
*            displays blocked, invariant, and       *
*                        variable-range reactions with their    *
*                        flux ranges                            *
***************************************************************** 
*To suppress output, output.lst=/dev/null
*To disable screen output w/o disabling .lst file creation, use command line w/ gams file.gms o /dev/null

***** Specifying the directories where the root files are present
$INLINECOM /*  */

***** Specifying the theoretical maximum biomass flux

$include "./min_flux_violation_GAMS_settings.txt"
* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal nscale 1e3

* Optional constraint on allowed proteome allocation to mitochondrial proteins (disable by setting to 1)
$setGlobal max_allowed_mito_proteome_allo_fraction 0.065

* Ribosome efficiency
$setGlobal kribonuc 13.2*3600
$setGlobal kribomito 13.2*3600

* Enforce part of proteome allocate to non-modeled protein
$setGlobal nonmodeled_proteome_allocation 0.45

options
    LP = soplex /*Solver selection*/
    limrow = 0 /*number of equations listed, 0 is suppresed*/
    limcol = 0 /*number of variables listed, 0 is suppresed*/
    iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
    decimals = 8 /*decimal places for display statement*/
    reslim = 1000000 /*wall-clock time limit for solver in seconds*/
    sysout = on /*solver status file report option*/
    solprint = on /*solution printing option*/
        
Sets
i
$include "%species_path%"
j
$include "%rxns_path%"
prosyn(j)
$include "%model_root_path%GAMS/model/RBA_rxns_prosyn.txt"
prowaste(j)
$include "%model_root_path%GAMS/model/RBA_rxns_prowaste.txt"
nuc_translation(j)
$include "%nuc_trans_path%"
mito_translation(j)
$include "%mito_trans_path%"
uptake(j) /*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"
media(j) /*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"
rxns_biomass(j)
$include "%biomass_path%"
*metabolism(j)
*$include "%metabolism_path%"
*enzsyn(j)
*$include "%enzsyn_path%"
*enzload(j)
*$include "%enzload_path%"
ribosyn(j)
$include "%model_root_path%GAMS/model/RBA_rxns_ribosyn.txt"

dummy(j)

;

****************************** PARAMETRIC CONSTANTS USED ********************************
*
*       s(i,j)          = Stoichiometry of metabolite i in reaction j
*       rxntype(j)      = Specified whether a reaction is irreversible (0), reversible
*                         (1 - forward, 2 - backward reaction), or exchange reactions (4)                  
*       LB(j)/UB(j)     = Stores the lower and upper bounds for each reaction j
*       low(j)/high(j)  = Stores the minimum/maximum flux od each reaction 
*       epsilon         = A small value
*       x               = counter
*
*****************************************************************************************

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
*kapp(j)
*$include "%kapp_path%"

LB(j) , UB(j)

high(j) , low(j)

x

epsilon /1e-6/
;

Variables

z
v(j)
;


*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e3 * %nscale%;

* Simulation top-level settings
* Enable or disable wasteful protein production, disabled by default (to solve faster)
* Note in solving: Enable protein waste flux might cause error for solver
* This is because enabling protein waste introduces several thousands more free variable to the system
* Thus, protein waste should only be implemented with actual data to constrain the free variable
*v.fx(j)$prowaste(j) = 0;

* Disable all biomass reactions
* Condition-specific biomass reaction activation has to be done in phenotype.txt file
v.fx(j)$rxns_biomass(j) = 0;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1e3 * %nscale%;

$include "%phenotype_path%"
* Set all organism-specific and/or condition-specific parameters in phenotype.txt (e.g., NGAM and biomass composition)

Equations

Obj
Stoic
Con_bio   
;

* DEFNINING THE CONSTRAINTS *************************************************************
Obj..           z =e= sum(j$dummy(j),v(j));
Stoic(i)..      sum(j, S(i,j) * v(j) )  =e=  0 ;

*****************************************************************************************

***** DEFINING THE UPPER AND LOWER BOUNDS FOR ALL THE REACTIONS *************************
scalar vmax /1000/;

*** Turn off reactions that are not active under aerobic conditions
*LB(j)$(offaeroglucose(j)) = 0;
*UB(j)$(offaeroglucose(j)) = 0;

*** SETTING LOWER AND UPPER BOUNDS FOR THE FLUXES (variable bounds)
*v.lo(j) = LB(j);
*v.up(j) = UB(j);
*****************************************************************************************

************************ DECLARING THE MODEL with constraints****************************
Model fva
/
Obj
Stoic
/;

fva.optfile = 1;

******SOLVING THE Linear Programming (LP) problem turning of one reaction at a time*******

*** loop through each reaction; dummy(j) has one reaction in its set at each iteration
*** Next, maximize/minimize their flux and store them in high(j)/low(j)
alias(j1,j)
*for(x = 1 to card(j),

*Updates output while solving, so errors won't get rid of output
file ff /%system.FN%.txt/;
put ff;
* Clear old file
put "rxn",system.tab,"upper-bound"/;
* Append old file
ff.ap=1;

loop(j1$(prosyn(j1)),
    dummy(j)$(ord(j) ne ord(j1)) = no;
    dummy(j)$(ord(j) eq ord(j1)) = yes;
    display j1;
    
    Solve fva using lp maximizing z;
    If(fva.modelstat eq 1,
        high(j1) = z.l;
	low(j1) = 0;
        );

    put ff;
    put j1.tl:0:50,system.tab,high(j1):0:8/;
    putclose ff
    );
