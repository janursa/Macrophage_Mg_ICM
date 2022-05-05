###--------- NF-kB pathway------------: 1 ###
""" inputs: pIKK
"""
NF_kB_main_block = """
            ##----variable definition and initial conditions----##
            IKB = 435;
            NFKB = 277;
            IKB_n = 25;
            NFKB_n = 931;
            IKB_NFKB_n = 29;
            IKB_NFKB = 68307;
            pIKB = 143;
            ## dummy variables ##
            IKB_prod = 1;
            ##------------------parameter definition------------##
            k11 = 10; # NFKB disassociation
            Kd11 = 10; # NFKB disassociation
            k12 = 0.5; # degradation of pIKB
            k13 = 0.5; # Nuclear translocation of NFKB and vice versa
            k14 = 0.5; # Nuclear translocation of NFKB and vice versa
            k15 = 10; # IkBa mRNA transcription
            Kd12 = 10; # IkBa mRNA transcription
            k16 = 0.5; # IkBa mRNA cytoplasm translocation and vice versa
            k17 = 0.5; # IkBa mRNA cytoplasm translocation and vice versa
            k18 = 0.5; # NFKB-IkBa bond formation
            k19 = 0.5; # NFKB-IkBa bond formation
            k20 = 0.5; # IkBa degredation
            ##---------------formulations---------------------##
            IKB_NFKB + pIKK => pIKB + NFKB + pIKK; k11*IKB_NFKB*pIKK/(pIKK + Kd11); ## NFkB disassosciation in cytoplasm ##
            pIKB => deg; k12*pIKB; ## pIKB degradation ##
            NFKB => NFKB_n; k13*NFKB-k14*NFKB_n ; ## Nuclear translocation of NFKB and vice versa ##
            $IKB_prod + NFKB_n => IKB_n + NFKB_n; k15*NFKB_n/(NFKB_n+Kd12); ## IkBa mRNA transcription ##
            IKB -> IKB_n; k16*IKB - k17*IKB_n; ## IkBa mRNA cytoplasm translocation and vice versa ##
            NFKB + IKB -> IKB_NFKB; k18*NFKB*IKB - k19*IKB_NFKB; ## NFKB-IkBa bond formation ##
            IKB => deg; k20*IKB; ## IkBa degredation ##
        """
free_params = dict( 
    k11 = [0,1000000], # NFKB disassociation
    Kd11 = [0,1000000], # NFKB disassociation
    k12 = [0,1], # degradation of pIKB
    k13 = [0,1], # Nuclear translocation of NFKB and vice versa
    k14 = [0,1], # Nuclear translocation of NFKB and vice versa
    k15 = [0,1000000], # IkBa mRNA transcription
    Kd12 = [0,1000000], # IkBa mRNA transcription
    k16 = [0,1], # IkBa mRNA cytoplasm translocation and vice versa
    k17 = [0,1], # IkBa mRNA cytoplasm translocation and vice versa
    k18 = [0,1], # NFKB-IkBa bond formation
    k19 = [0,1], # NFKB-IkBa bond formation
    k20 = [0,1], # IkBa degredation
)