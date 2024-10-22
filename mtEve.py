from Bio import Phylo
from io import StringIO

newick_tree = """((((((((((((((((((((((KY077676.1_Homo_sapiens_haplogroup_H7j1_mitochondrion_complete_genome:0.00030189,KY934476.1_Homo_sapiens_haplogroup_H1h1_mitochondrion_complete_genome:0.00030189)'-':0.00015110,KP340170.1_Homo_sapiens_isolate_azerE43_haplogroup_HV4a1_mitochondrion_complete_genome:0.00045300)'-':0.00011086,(JF343123.1_Homo_sapiens_haplogroup_V3b1_mitochondrion_complete_genome:0.00036233,KY496869.1_Homo_sapiens_isolate_MP_13_haplogroup_V_mitochondrion_complete_genome:0.00036233)'-':0.00020152)'-':0.00030064,JN084084.1_Homo_sapiens_haplogroup_N9a10_mitochondrion_complete_genome:0.00086450)'-':0.00005747,(KU521491.1_Homo_sapiens_isolate_KK13_haplogroup_Y2a1_mitochondrion_complete_genome:0.00006035,KU521494.1_Homo_sapiens_isolate_MND64_haplogroup_Y2a1_mitochondrion_complete_genome:0.00006035)'-':0.00086161)'-':0.00002664,KT819263.1_Homo_sapiens_isolate_Fig18_haplogroup_L3e5a1_mitochondrion_complete_genome:0.00094860)'-':0.00004074,(HM448049.1_Homo_sapiens_haplogroup_X2d_mitochondrion_complete_genome:0.00054373,HM453712.1_Homo_sapiens_haplogroup_X2B_mitochondrion_complete_genome:0.00054373)'-':0.00044561)'-':0.00003060,KU521454.1_Homo_sapiens_isolate_PAI17_haplogroup_F3b1_mitochondrion_complete_genome:0.00101994)'-':0.00001578,HM804485.1_Homo_sapiens_haplogroup_U6a7a1_mitochondrion_complete_genome:0.00103571)'-':0.00000638,KY303770.1_Homo_sapiens_haplogroup_R2c_mitochondrion_complete_genome:0.00104209)'-':0.00005974,(KX459697.1_Homo_sapiens_haplogroup_A4-A200G_mitochondrion_complete_genome:0.00066481,JQ247408.1_Homo_sapiens_haplogroup_A2f1a_mitochondrion_complete_genome:0.00066481)'-':0.00043703)'-':0.00002869,KY686210.1_Homo_sapiens_haplogroup_R6_mitochondrion_complete_genome:0.00113053)'-':0.00001789,(HQ189135.1_Homo_sapiens_haplogroup_K2a3_mitochondrion_complete_genome:0.00057398,KT698008.1_Homo_sapiens_isolate_57_Sb_haplogroup_K1a4a1_mitochondrion_complete_genome:0.00057398)'-':0.00057444)'-':0.00004242,(KY934478.1_Homo_sapiens_haplogroup_W6_mitochondrion_complete_genome:0.00036235,KU508374.1_Homo_sapiens_haplogroup_W1c_mitochondrion_complete_genome:0.00036235)'-':0.00082849)'-':0.00002924,(KM986616.1_Homo_sapiens_isolate_Y488_haplogroup_N1a1a3_mitochondrion_complete_genome:0.00087657,(KY348642.1_Homo_sapiens_haplogroup_I1a1a3_mitochondrion_complete_genome:0.00058900,(KY369152.2_Homo_sapiens_haplogroup_I4a_mitochondrion_complete_genome:0.00039251,GU590993.1_Homo_sapiens_haplogroup_I3_mitochondrion_complete_genome:0.00039251)'-':0.00019649)'-':0.00028757)'-':0.00034350)'-':0.00004520,(((JN419195.1_Homo_sapiens_haplogroup_T2b16_mitochondrion_complete_genome:0.00030194,KX440315.1_Homo_sapiens_isolate_JT129_haplogroup_T2a1a_mitochondrion_complete_genome:0.00030194)'-':0.00018136,JF837819.1_Homo_sapiens_haplogroup_T1a1_mitochondrion_complete_genome:0.00048330)'-':0.00068282,(HQ914447.1_Homo_sapiens_haplogroup_J1b_mitochondrion_complete_genome:0.00083127,(KX440275.1_Homo_sapiens_isolate_JT89_haplogroup_J2b1_mitochondrion_complete_genome:0.00069503,KX440262.1_Homo_sapiens_isolate_JT76_haplogroup_J2a2d1_mitochondrion_complete_genome:0.00069503)'-':0.00013624)'-':0.00033485)'-':0.00009915)'-':0.00011042,(KM986608.1_Homo_sapiens_isolate_Y453_haplogroup_L4b2a1_mitochondrion_complete_genome:0.00133126,(KY686209.1_Homo_sapiens_haplogroup_M60_mitochondrion_complete_genome:0.00108872,KM986625.1_Homo_sapiens_isolate_Y539_haplogroup_M23_mitochondrion_complete_genome:0.00108872)'-':0.00024254)'-':0.00004444)'-':0.00007449,KT698006.1_Homo_sapiens_isolate_52_Sb_haplogroup_U2e2a4_mitochondrion_complete_genome:0.00145018)'-':0.00051737,(KR135861.1_Homo_sapiens_isolate_SUD103_haplogroup_L2e1_mitochondrion_complete_genome:0.00154372,KR135883.1_Homo_sapiens_isolate_MOZ326_haplogroup_L2a1a_mitochondrion_complete_genome:0.00154372)'-':0.00042384)'-':0.00061419,FJ713601.1_Homo_sapiens_haplogroup_L1c1d_mitochondrion_complete_genome:0.00258174)'-':0.00380994,((KX198087.1_Homo_sapiens_neanderthalensis_isolate_GoyetQ305-4_mitochondrion_complete_genome:0.00033227,KX198088.1_Homo_sapiens_neanderthalensis_isolate_GoyetQ57-2_mitochondrion_complete_genome:0.00033227)'-':0.00007566,(KX198085.1_Homo_sapiens_neanderthalensis_isolate_GoyetQ374a-1_mitochondrion_complete_genome:0.00000000,(KX198086.1_Homo_sapiens_neanderthalensis_isolate_GoyetQ305-7_mitochondrion_complete_genome:0.00000000,KX198084.1_Homo_sapiens_neanderthalensis_isolate_GoyetQ56-1_mitochondrion_complete_genome:0.00000000)'-':0.00000000)'-':0.00040793)'-':0.00598375)'-':0.00000000,(KT780370.1_Homo_sapiens_ssp._Denisova_isolate_Denisova8_mitochondrion_complete_genome:0.00260896,(FR695060.1_Homo_sp._Altai_complete_mitochondrial_genome_isolate_Denisova_molar:0.00006035,FN673705.1_Homo_sp._Altai_complete_mitochondrial_genome_sequence_from_Denisova_Altai_Russia:0.00006035)'-':0.00254861)'-':0.01443023);
"""
tree = Phylo.read(StringIO(newick_tree), "newick")

mutation_rate_lower = 1.26e-8  
mutation_rate_upper = 2.8e-8   

#Branch length between two specific clades and their MRCA
def separation_time(clade1, clade2, mutation_rate):
    mrca = tree.common_ancestor(clade1, clade2)
    branch_length_clade1 = mrca.distance(clade1)
    branch_length_clade2 = mrca.distance(clade2)   
    total_branch_length = (branch_length_clade1 + branch_length_clade2)/2  #2 branch lenghts used
    separation_time = total_branch_length / mutation_rate
    return separation_time

#Example taxa names for Neanderthals and Denisova
neanderthal_clades = [
    "KX198087.1_Homo_sapiens_neanderthalensis_isolate_GoyetQ305-4_mitochondrion_complete_genome"
]

denisovan_clades = [
    "KT780370.1_Homo_sapiens_ssp._Denisova_isolate_Denisova8_mitochondrion_complete_genome"
]

#Separation time for Neanderthals and Denisova using their respective clades
separation_time_lower = int(separation_time(neanderthal_clades[0], denisovan_clades[0], mutation_rate_lower))
separation_time_upper = int(separation_time(neanderthal_clades[0], denisovan_clades[0], mutation_rate_upper))
separation_time_average = int(separation_time_upper + (separation_time_lower - separation_time_upper)/2)

print(f"Estimated separation time between Neanderthals and Denisova (lower bound): {separation_time_lower} years ago")
print(f"Estimated separation time between Neanderthals and Denisova (upper bound): {separation_time_upper} years ago")
print(f"Estimated separation time between Neanderthals and Denisova (average bound): {separation_time_average} years ago")
#mtEve age calculation
def calculate_mitochondrial_eve_age(tree, mutation_rate):
    root = tree.root
    total_branch_length = sum(clade.branch_length for clade in root.get_terminals())/12 #12 branch lenghts which > 0.001. Others are significantly lower and hhave no weight in calculations
    eve_age = total_branch_length / mutation_rate
    return eve_age

mitochondrial_eve_age_lower = int(calculate_mitochondrial_eve_age(tree, mutation_rate_lower))
mitochondrial_eve_age_upper = int(calculate_mitochondrial_eve_age(tree, mutation_rate_upper))
mitochondrial_eve_age_average = int(mitochondrial_eve_age_upper + (mitochondrial_eve_age_lower - mitochondrial_eve_age_upper)/2)


print(f"Mitochondrial Eve age (lower bound): {mitochondrial_eve_age_lower} years ago")
print(f"Mitochondrial Eve age (upper bound): {mitochondrial_eve_age_upper} years ago")
print(f"Mitochondrial Eve age (average bound): {mitochondrial_eve_age_average} years ago")
