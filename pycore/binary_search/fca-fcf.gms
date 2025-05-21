***************************************************************************************
*                                CHE 597C- HW #3                                      *
*                  FBA for the iAF1260 metabolic model of E. coli                     *
***************************************************************************************

$INLINECOM /*  */
$onlisting
$offdigit

$set biomass Ec_biomass_iAF1260_core_59p81M
$setGlobal metabolism_path ../../GAMS/model/RBA_rxns_rxnmetabolicnetwork.txt
$setGlobal metabolism_path ../../GAMS/model/RBA_rxns_rxnmetabolicnetwork.txt

OPTIONS 
       decimals = 8
       solprint = on
       reslim = 1000000
       iterlim = 10000000
       domlim = 10
       limcol = 1000
       limrow = 1000
       optca = 0.0
       optcr = 1E-9
       mip = cplex
;

SETS
        i Set of metabolites 
$include metabolites.txt 

        j Set of reactions 
$include reactions.txt 

        reversible(j) Reversible rxns (other than exchange rxns)
$include reversible_reactions_no_exchange.txt 

        irreversible(j)  Irreversible rxns
$include irreversible_reactions.txt 

        exchange(j) Exchange rxns
$include exchange_reactions.txt 

        medium(j) Exchange rxns corresponding to compounds present in the growth medium 
$include minimal_medium.txt 

        regulation(j) Reactions that should off due to regulatory constraints 
$include regulated_reactions.txt 

        blocked(j) The set of blocked rxns
$include Blocked_found.txt

        current_num(j) The current flux in numerator
 
        current_den(j) The current flux in denominator 

        coupled(j) The set of reactions coupled with j
;

PARAMETERS
  UB(j) Lowerbound on reaction fluxes 

  LB(j) Upperbound on reaction fluxes 

  S(i,j) contains the Stoichiometric matrix of the metabolic model
$include S_matrix.txt 

  p(j)  A parameter for biomass

  Rmin

  Rmax

;

p('Ec_biomass_iAF1260_core_59p81M') = 1;

***** Set the bounds *****
LB(j)$irreversible(j) = 0;
UB(j)$irreversible(j) = 1000;

* Reversible reactions 
LB(j)$reversible(j) = 0;
UB(j)$reversible(j) = 1000;

* Exchange reactions
LB(j)$exchange(j) = 0;
UB(j)$exchange(j) = 1000;

***** Set the conditions for the growth medium *****
UB(j)$medium(j) = 1000;

UB('EX_glc(e)') = 1000;
UB('EX_o2(e)') = 1000;

* Turn off all reactions in the set regulation
UB(j)$(regulation(j)) = 0;

* Turn off the wild-type biomass equation
UB('Ec_biomass_iAF1260_WT_59p81M') = 0;

POSITIVE VARIABLES
        v(j)      Flux 
        r         Normalization factor
;

*v.up(j)=UB(j);
*v.lo(j)=LB(j);

VARIABLES
        z         Objective function 
;

EQUATIONS
        obj              Objective function 
        massbalance(i)   Mass balance equations for each metabolite i
        den_const        The flux in denominator of the FCF
        glc_const        Constraints on glucose uptake
        o2_const         Constraint on oxygen uptake
        ATPM_const       Constraint on ATPM
        con1a(j)
        con1b(j)
;

obj..             z =e= v('%biomass%');
massbalance(i)..  sum(j,S(i,j)*v(j)) =e= 0;
den_const..       sum(j$current_den(j),v(j)) =e= 1;
*den_const..       v(j)$current_den(j) =e= 1;
*glc_const..       v('EX_glc(e)') =g= -10*r;
*o2_const..        v('EX_o2(e)') =g= -20*r;
*ATPM_const..      v('ATPM') =e= 8.39*r;
con1a(j)..        v(j)  =g=  LB(j)*r;
con1b(j)..      v(j)  =l=  UB(j)*r;

************** Model definitions ********************
Model FCF 
/
  obj
  massbalance
  den_const
*  glc_const
*  o2_const
*  ATPM_const
/;

FCF.optfile = 1;


FILE res /coupling_results_biomass.txt/;
PUT res;

alias(j,j1,j2);

LOOP(j2$(not blocked(j2)),
    current_den(j) = no;
    current_den(j2) = yes;

    SOLVE FCF USING LP MAXIMIZING z;
    if(FCF.modelstat = 1,
      Rmax = z.l;
    elseif (FCF.modelstat = 3),
      Rmax = -1;    
    );

    if([FCF.modelstat = 1] or [FCF.modelstat = 3],
       SOLVE FCF USING LP MINIMIZING z;
       if(FCF.modelstat = 1,
          Rmin = z.l;

          if([Rmin = 0] and [Rmax > 0],
             PUT "biomass --> ",j2.tl:0,"  (Rmin , Rmax) = (",Rmin:0:8,",",Rmax:0:8,")"/;
          elseif ([Rmin > 0] and [Rmax > 0] and [Rmax > Rmin]), 
             PUT "biomass <--> ",j2.tl:0,"  (Rmin , Rmax) = (",Rmin:0:8,",",Rmax:0:8,")"/;
          elseif ([Rmin > 0] and [Rmax > 0] and [Rmin = Rmax]), 
             PUT "biomass <==> ",j2.tl:0,"  (Rmin , Rmax) = (",Rmin:0:8,",",Rmax:0:8,")"/;
          elseif ([Rmin > 0] and [Rmax < 0]), 
             PUT j2.tl:0," --> biomass  (Rmin , Rmax) = (",Rmin:0:8,",",Rmax:0:8,")"/;
          );

       else
         PUT "Error2! ",j2.tl:0," model status = ",FCF.modelstat:0:0/;
       );
   else
       PUT "Error1! ",j2.tl:0," model status = ",FCF.modelstat:0:0/;
   );
);


