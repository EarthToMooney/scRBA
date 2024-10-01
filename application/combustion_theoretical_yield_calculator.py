# combusting 1 mol of product to CO2 and H2O, and doing the same for the reactant, to determine yields
# helps test in a model-agnostic way if yields makes sense
# these estimates are maximum yields, assuming you're not taking up carbons from elsewhere
import chempy
def parse_formula_charge(metabolite):
    try:
        # here we take the charged form
        charge = metabolite['charges']
        if len(charge) > 1:
            #print ('multi')
            if charge[0]:
                charge = charge[0]
                formula = metabolite['formulae'][0]
            else:
                charge = charge[1]
                if len(metabolite['formulae']) > 1:
                    formula = metabolite['formulae'][1]
                else:
                    formula = metabolite['formulae'][0]
        else:
            charge = metabolite['charges'][0]
            formula = metabolite['formulae'][0]
    except:
        print ("Error with passed metabolite")
        print (metabolite)#print ("\t",formula, charge)
    return formula, charge
 
def combustion(metabolite):
    try:
        # here we take the charged form
        formula, charge = parse_formula_charge(metabolite)
        if charge:
            a = chempy.Substance.from_formula(str(formula)+str(charge))
        else:
            a = chempy.Substance.from_formula(str(formula))
        # compute free electrons by balancing complete combusion of CxHyOz
        x = a.composition[6] #C, carbon
        y = a.composition[1] #H, hydrogen
        z = a.composition[8] #O, oxygen
        e = 4*x + y - 2*z - charge
        #print (a)
    except:
        print ("Error with passed metabolite")
        print (metabolite)
    return e
 
met = '3hpp'
sub = bigg_metabolites[met]['formulae'][0]+str(bigg_metabolites[met]['charges'][0])
print(sub)
b = chempy.Substance.from_formula(sub)

reac, prod = chempy.balance_stoichiometry({bigg_metabolites['succ']['formulae'][0], 'O'}, {'H2O', 'CO2'})
print(reac, prod)

def balance(r,p):
    try:
        reac,prod = chempy.balance_stoichiometry({parse_formula_charge(r)[0], 'O'},
                                                 {parse_formula_charge(p)[0],'H2O'},
                                                 underdetermined=None)
        status = 0
    except ValueError:
        print('Error with',acid,item)
        reac,prod = chempy.balance_stoichiometry({parse_formula_charge(r)[0], 'O','H'},
                                                 {parse_formula_charge(p)[0],'H2O'},
                                                 underdetermined=None)
        #reac,prod = chempy.balance_stoichiometry({parse_formula_charge(r)[0], 'O','HCO3'},
        #                                  {parse_formula_charge(p)[0],'H2O','H','CO2'},
        #                                 underdetermined=None)
        status = 1
        #reac,prod = {parse_formula_charge(r)[0]:1},{parse_formula_charge(p)[0]:0}
    return(reac,prod,status)