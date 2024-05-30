class GAMS_settings:
    def __init__(self,
                 species_path='../model/RBA_species.txt',
                 rxns_path='../model/RBA_rxns.txt',
                 prosyn_path='../model/RBA_rxns_prosyn.txt',
                 nuc_trans_path='../model/RBA_nuc_translation.txt',
                 mito_trans_path='../model/RBA_mito_translation.txt',
                 uptake_path='../model/RBA_rxns_EXREV.txt',
                 media_path='../model/RBA_rxns_EXREV_YNB.txt',
                 sij_path='../model/RBA_sij.txt',
                 prolen_path='../model/RBA_proteinLength.txt',
                 kapp_path='../model/RBA_kapp.txt',
                 enz_cap_declares_path='../model/RBA_enzCapacityConstraints_declares.txt',
                 enz_cap_eqns_path='../model/RBA_enzCapacityConstraints_eqns.txt',
                 kribonuc='10.5*3600',
                 kribomito='10.5*3600'):
        self.species_path = species_path
        self.rxns_path = rxns_path
        self.prosyn_path = prosyn_path
        self.nuc_trans_path = nuc_trans_path
        self.mito_trans_path = mito_trans_path
        self.uptake_path = uptake_path
        self.media_path = media_path
        self.sij_path = sij_path
        self.prolen_path = prolen_path
        self.kapp_path = kapp_path
        self.enz_cap_declares_path = enz_cap_declares_path
        self.enz_cap_eqns_path = enz_cap_eqns_path
        self.kribonuc = kribonuc
        self.kribomito = kribomito
    
    def export_to_txt_file(self, filepath):
        props = ['species_path', 'rxns_path', 'prosyn_path', 'nuc_trans_path',
                 'mito_trans_path', 'uptake_path', 'media_path', 'sij_path',
                 'prolen_path', 'kapp_path', 'enz_cap_declares_path',
                 'enz_cap_eqns_path', 'kribonuc', 'kribomito']
        text = []
        for p in props:
            text.append('$setGlobal ' + p + ' ' + self.__dict__[p])
        
        with open(filepath, 'w') as f:
            f.write('\n'.join(text))
            
    def import_from_txt_file(self, filepath):
        with open(filepath) as f:
            text = f.read().split('\n')
        text = [i[11:] for i in text if i != '']
        for i in text:
            k,v = i.split(' ')
            self.__setattr__(k,v)
