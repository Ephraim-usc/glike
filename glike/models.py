import msprime
from .glike import *


########## Twoway Admixture ###########
def twoway_admixture_demo(t1, t2, t3, r1, N, N_a, N_b, N_b_, N_e):
  demo = Demo()
  demo.add_phase(Phase(0, t1, [1/N]))
  demo.add_phase(Phase(t1, t2, [1/N_a, 1/N_b], P = np.array([[r1, 1-r1]])))
  demo.add_phase(Phase(t2, t3, [1/N_a, 1/N_b_]))
  demo.add_phase(Phase(t3, np.inf, [1/N_e], P = np.array([[1], [1]])))
  return demo

def twoway_admixture_demography(t1, t2, t3, r1, N, N_a, N_b, N_b_, N_e):
  demography = msprime.Demography()
  demography.add_population(name = "O", initial_size = N)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.add_population(name = "E", initial_size = N_e)
  demography.add_admixture(time=t1, derived="O", ancestral=["A", "B"], proportions = [r1, 1-r1])
  demography.add_population_parameters_change(time = t2, initial_size = N_b_, population = "B")
  demography.add_population_split(time=t3, derived=["A", "B"], ancestral="E")
  
  return demography


########## Threeway Split ###########
def threeway_split_demo(t1, t2, N_a, N_b, N_c, N_d, N_e):
  demo = Demo()
  demo.add_phase(Phase(0, t1, [1/N_a, 1/N_b, 1/N_c], populations = ["A", "B", "C"]))
  demo.add_phase(Phase(t1, t2, [1/N_d, 1/N_c], P = np.array([[1, 0], [1, 0], [0, 1]])))
  demo.add_phase(Phase(t2, math.inf, [1/N_e], P = np.array([[1], [1]])))
  return demo

def threeway_split_demography(t1, t2, N_a, N_b, N_c, N_d, N_e):
  demography = msprime.Demography()
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.add_population(name = "C", initial_size = N_c)
  demography.add_population(name = "D", initial_size = N_d)
  demography.add_population(name = "E", initial_size = N_e)
  
  demography.add_population_split(time=t1, derived=["A", "B"], ancestral="D")
  demography.add_population_split(time=t2, derived=["D", "C"], ancestral="E")
  
  return demography


########## Threeway Admixture ###########
def threeway_admixture_demo(t1 = 30, t2 = 60, t3 = 10000, r1 = 0.4, r2 = 0.7, N = 2000, N_a = 20000, N_b = 3000, N_c = 30000, N_d = 10000, N_e = 5000):
  demo = Demo()
  demo.add_phase(Phase(0, t1, [1/N]))
  demo.add_phase(Phase(t1, t2, [1/N_a, 1/N_b], P = np.array([[r1, 1-r1]])))
  demo.add_phase(Phase(t2, t3, [1/N_a, 1/N_c, 1/N_d], P = np.array([[1, 0, 0], [0, r2, 1-r2]])))
  demo.add_phase(Phase(t3, math.inf, [1/N_e], P = np.array([[1], [1], [1]])))
  return demo

def threeway_admixture_demography(t1 = 30, t2 = 60, t3 = 10000, r1 = 0.4, r2 = 0.7, N = 2000, N_a = 20000, N_b = 3000, N_c = 30000, N_d = 10000, N_e = 5000):
  demography = msprime.Demography()
  demography.add_population(name = "O", initial_size = N)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.add_population(name = "C", initial_size = N_c)
  demography.add_population(name = "D", initial_size = N_d)
  demography.add_population(name = "E", initial_size = N_e)
  
  demography.add_admixture(time=t1, derived="O", ancestral=["A", "B"], proportions = [r1, 1-r1])
  demography.add_admixture(time=t2, derived="B", ancestral=["C", "D"], proportions = [r2, 1-r2])
  demography.add_population_split(time=t3, derived=["A", "C", "D"], ancestral="E")
  
  return demography


########## Neandertal admixture (stdpopsim 3I21) ###########
def neandertal_admixture_demo(t1, t2, t3, t4, N_yri, N_ceu, N_nea, m1):
  Q = np.array([[0, 0, 0], [0, 0, 0], [0, m1, -m1]])
  demo = Demo()
  demo.add_phase(Phase(0, t1, [1/N_yri, 1/N_ceu, 1/N_nea], populations = ["yri", "ceu", "nea"]))
  demo.add_phase(Phase(t1, t2, [1/N_yri, 1/N_ceu, 1/N_nea], Q = Q, populations = ["yri", "ceu", "nea"]))
  demo.add_phase(Phase(t2, t3, [1/N_yri, 1/N_ceu, 1/N_nea], populations = ["yri", "ceu", "nea"]))
  demo.add_phase(Phase(t3, t4, [1/N_yri, 1/N_nea], P = np.array([[1, 0], [1, 0], [0, 1]]), populations = ["yri", "nea"]))
  demo.add_phase(Phase(t4, math.inf, [1/N_yri], P = np.array([[1], [1]]), populations = ["yri"]))
  return demo

def neandertal_admixture_demography(t1, t2, t3, t4, N_yri, N_ceu, N_nea, m1):
  demography = msprime.Demography()
  demography.add_population(name = "yri", initial_size = N_yri)
  demography.add_population(name = "ceu", initial_size = N_ceu)
  demography.add_population(name = "nea", initial_size = N_nea)
  demography.add_migration_rate_change(time = t1, rate = m1, source = "nea", dest = "ceu")
  demography.add_migration_rate_change(time = t2, rate = 0, source = "nea", dest = "ceu")
  demography.add_mass_migration(t3, source = "ceu", dest = "yri", proportion = 1)
  demography.add_mass_migration(t4, source = "nea", dest = "yri", proportion = 1)
  return demography



########## American Admixture (stdpopsim 4B11) ###########
def american_admixture_demo(t1 = 12, t2 = 920, t3 = 2040, t4 = 5920, 
                            r1 = 0.17, r2 = 0.33, 
                            N_afr = 14474, N_eur = 34038, N_asia = 45851, N_admix = 54663, N_ooa = 1867, N_anc = 7310,
                            gr_eur = 0.0038, gr_asia = 0.0048, gr_admix = 0.05,
                            m1 = 0.0, m2 = 0.0, m3 = 0.0, m4 = 0.0):
  Q0 = np.array([[-m1-m2, m1, m2, 0], [m1, -m1-m3, m3, 0], [m2, m3, -m2-m3, 0], [0, 0, 0, 0]]) if max(m1, m2, m3) > 0 else None
  Q1 = np.array([[-m1-m2, m1, m2], [m1, -m1-m3, m3], [m2, m3, -m2-m3]]) if max(m1, m2, m3) > 0 else None
  Q2 = np.array([[-m4, m4], [m4, -m4]]) if m4 > 0 else None
  
  demo = Demo()
  demo.add_phase(Phase(0, 1e-6, [1/N_afr, 1/N_eur, 1/N_asia, 1/N_admix], grs = [0, gr_eur, gr_asia, gr_admix], populations = ["afr", "eur", "asia", "admix"]))
  demo.add_phase(Phase(1e-6, t1, [1/N_afr, 1/N_eur, 1/N_asia, 1/N_admix], grs = [0, gr_eur, gr_asia, gr_admix], Q = Q0, populations = ["afr", "eur", "asia", "admix"]), discretize = 200)
  
  P_admixture = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [r1, r2, 1-r1-r2]])
  demo.add_phase(Phase(t1, t2, [1/N_afr, 1/N_eur*math.exp(gr_eur * t1), 1/N_asia*math.exp(gr_asia * t1)], [0, gr_eur, gr_asia], P = P_admixture, Q = Q1, populations = ["afr", "eur", "asia"]), discretize = 200)
  
  P_asia_split = np.array([[1, 0], [0, 1], [0, 1]])
  demo.add_phase(Phase(t2, t3, [1/N_afr, 1/N_ooa], P = P_asia_split, Q = Q2, populations = ["afr", "eur"]), discretize = 500)
  
  P_ooa_split = np.array([[1], [1]])
  demo.add_phase(Phase(t3, t4, [1/N_afr], P = P_ooa_split, populations = ["afr"]))
  
  demo.add_phase(Phase(t4, math.inf, [1/N_anc], populations = ["afr"]))
  return demo

def american_admixture_demography(t1 = 12, t2 = 920, t3 = 2040, t4 = 5920, 
                                  r1 = 0.17, r2 = 0.33, 
                                  N_afr = 14474, N_eur = 34038, N_asia = 45851, N_admix = 54663, N_ooa = 1867, N_anc = 7310,
                                  gr_eur = 0.0038, gr_asia = 0.0048, gr_admix = 0.05,
                                  m1 = 0.0, m2 = 0.0, m3 = 0.0, m4 = 0.0):
  demography = msprime.Demography()
  demography.add_population(name = "afr", initial_size = N_afr)
  demography.add_population(name = "eur", initial_size = N_eur, growth_rate = gr_eur)
  demography.add_population(name = "asia", initial_size = N_asia, growth_rate = gr_asia)
  demography.add_population(name = "admix", initial_size = N_admix, growth_rate = gr_admix)
  demography.set_symmetric_migration_rate(populations = ["afr", "eur"], rate = m1)
  demography.set_symmetric_migration_rate(populations = ["afr", "asia"], rate = m2)
  demography.set_symmetric_migration_rate(populations = ["eur", "asia"], rate = m3)
  
  demography.add_admixture(time=t1, derived="admix", ancestral=["afr", "eur", "asia"], proportions = [r1, r2, 1-r1-r2])
  
  demography.add_mass_migration(time=t2, source="asia", dest="eur", proportion=1)
  demography.add_population_parameters_change(time=t2, initial_size = N_ooa, growth_rate=0, population="eur")
  demography.add_symmetric_migration_rate_change(time = t2, populations = ["afr", "eur"], rate = m4)
  
  demography.add_mass_migration(time=t3, source="eur", dest="afr", proportion=1)
  demography.add_symmetric_migration_rate_change(time = t3, populations = ["afr", "eur"], rate = 0)
  
  demography.add_population_parameters_change(time=t4, initial_size = N_anc, growth_rate=0, population="afr")
  return demography


########## Ancient Europe (stdpopsim 4A21) ###########
def ancient_europe_demo(t1, t2, t3, t4, t5, t6, r1, r2, r3, N_ana, N_neo, N_whg, N_bronze, N_yam, N_ehg, N_chg, N_ne, N_wa, N_ooa, gr):
  demo = Demo()
  demo.add_phase(Phase(0, t1, 
                       ns = [1/N_ana, 1/N_neo, 1/N_whg, 1/(N_bronze*math.exp(gr*t1)), 1/N_yam, 1/N_ehg, 1/N_chg],
                       grs = [0.0, 0.0, 0.0, gr, 0.0, 0.0, 0.0],
                       populations = ["ana", "neo", "whg", "bronze", "yam", "ehg", "chg"]))
  P_bronze_admixture = np.array([
      [1, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [0, 0, 1, 0, 0, 0],
      [0, r1, 0, 1-r1, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1]
  ])
  demo.add_phase(Phase(t1, t2, [1/N_ana, 1/N_neo, 1/N_whg, 1/N_yam, 1/N_ehg, 1/N_chg], P = P_bronze_admixture,
                      populations = ["ana", "neo", "whg", "yam", "ehg", "chg"]))
  P_yam_admixture = np.array([
      [1, 0, 0, 0, 0],
      [0, 1, 0, 0, 0],
      [0, 0, 1, 0, 0],
      [0, 0, 0, r2, 1-r2],
      [0, 0, 0, 1, 0],
      [0, 0, 0, 0, 1]
  ])
  demo.add_phase(Phase(t2, t3, [1/N_ana, 1/N_neo, 1/N_whg, 1/N_ehg, 1/N_chg], P = P_yam_admixture,
                      populations = ["ana", "neo", "whg", "ehg", "chg"]))
  P_neo_admixture = np.array([
      [1, 0, 0, 0],
      [r3, 1-r3, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1]
  ])
  demo.add_phase(Phase(t3, t4, [1/N_ana, 1/N_whg, 1/N_ehg, 1/N_chg], P = P_neo_admixture,
                      populations = ["ana", "whg", "ehg", "chg"]))
  P_ne_split = np.array([
      [1, 0, 0],
      [0, 1, 0],
      [0, 1, 0],
      [0, 0, 1]
  ])
  demo.add_phase(Phase(t4, t5, [1/N_ana, 1/N_ne, 1/N_chg], P = P_ne_split,
                      populations = ["ana", "ne", "chg"]))
  P_wa_split = np.array([
      [1, 0],
      [0, 1],
      [1, 0]
  ])
  demo.add_phase(Phase(t5, t6, [1/N_wa, 1/N_ne], P = P_wa_split,
                      populations = ["wa", "ne"]))
  P_ooa_split = np.array([
      [1],
      [1]
  ])
  demo.add_phase(Phase(t6, math.inf, [1/N_ooa], P = P_ooa_split,
                      populations = ["ooa"]))
  return demo

def ancient_europe_demography(t1, t2, t3, t4, t5, t6, r1, r2, r3, N_ana, N_neo, N_whg, N_bronze, N_yam, N_ehg, N_chg, N_ne, N_wa, N_ooa, gr):
  demography = msprime.Demography()
  demography.add_population(name = "ana", initial_size = N_ana)
  demography.add_population(name = "neo", initial_size = N_neo)
  demography.add_population(name = "whg", initial_size = N_whg)
  demography.add_population(name = "bronze", initial_size = N_bronze, growth_rate = gr)
  demography.add_population(name = "yam", initial_size = N_yam)
  demography.add_population(name = "ehg", initial_size = N_ehg)
  demography.add_population(name = "chg", initial_size = N_chg)
  demography.add_population(name = "ne", initial_size = N_ne)
  demography.add_population(name = "wa", initial_size = N_wa)
  demography.add_population(name = "ooa", initial_size = N_ooa)
  
  demography.add_admixture(time=t1, derived="bronze", ancestral=["neo", "yam"], proportions = [r1, 1-r1])
  demography.add_admixture(time=t2, derived="yam", ancestral=["ehg", "chg"], proportions = [r2, 1-r2])
  demography.add_admixture(time=t3, derived="neo", ancestral=["ana", "whg"], proportions = [r3, 1-r3])
  
  demography.add_population_split(time=t4, derived=["whg", "ehg"], ancestral="ne")
  demography.add_population_split(time=t5, derived=["ana", "chg"], ancestral="wa")
  demography.add_population_split(time=t6, derived=["ne", "wa"], ancestral="ooa")
  return demography


########## Native Hawaiian ###########
def native_hawaiians_demo(t1 = 19, t2 = 411, t3 = 1040, t4 = 2004, r1 = 0.0, r2 = 0.198, r3 = 0.334, 
                          N_admix = 35682, N_afr = 10000, N_eur = 13388, N_asia = 25234, N_pol = 15695, N_aa = 2702, N_ooa = 2470, N_anc = 2665, gr = 0.078, 
                          m_afr_eur = 0.0, m_afr_asia = 0.0, m_afr_pol = 0.0, m_eur_asia = 0.0, m_eur_pol = 0.0, m_asia_pol = 0.0):
  Q = np.array([[-m_afr_eur-m_afr_asia-m_afr_pol, m_afr_eur, m_afr_asia, m_afr_pol], 
                [m_afr_eur, -m_afr_eur-m_eur_asia-m_eur_pol, m_eur_asia, m_eur_pol], 
                [m_afr_asia, m_eur_asia, -m_afr_asia-m_eur_asia-m_asia_pol, m_asia_pol],
                [m_afr_pol, m_eur_pol, m_asia_pol, -m_afr_pol-m_eur_pol-m_asia_pol]]) if max(m_afr_eur, m_afr_asia, m_afr_pol, m_eur_asia, m_eur_pol, m_asia_pol) > 0 else None
  
  demo = Demo()
  demo.add_phase(Phase(0, t1, [1/N_admix * math.exp(-t1 * gr)], [gr], populations = ["admix"]))
  P_admixture = np.array([
      [r1, r2, r3, 1-r1-r2-r3]
  ])
  demo.add_phase(Phase(t1, t2, [1/N_afr, 1/N_eur, 1/N_asia, 1/N_pol], P = P_admixture, Q = Q, populations = ["afr", "eur", "asia", "pol"]), discretize = 100)
  P_pol_split = np.array([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
      [0, 0, 1]
  ])
  demo.add_phase(Phase(t2, t3, [1/N_afr, 1/N_eur, 1/N_aa], P = P_pol_split, populations = ["afr", "eur", "asia"]))
  P_asia_split = np.array([
      [1, 0],
      [0, 1],
      [0, 1]
  ])
  demo.add_phase(Phase(t3, t4, [1/N_afr, 1/N_ooa], P = P_asia_split, populations = ["afr", "eur"]))
  P_ooa_split = np.array([
      [1],
      [1]
  ])
  demo.add_phase(Phase(t4, math.inf, [1/N_anc], P = P_ooa_split))
  return demo

def native_hawaiians_demography(t1 = 19, t2 = 411, t3 = 1040, t4 = 2004, r1 = 0.0, r2 = 0.198, r3 = 0.334, 
                                N_admix = 35682, N_afr = 10000, N_eur = 13388, N_asia = 25234, N_pol = 15695, N_aa = 2702, N_ooa = 2470, N_anc = 2665, gr = 0.078, 
                                m_afr_eur = 0.0, m_afr_asia = 0.0, m_afr_pol = 0.0, m_eur_asia = 0.0, m_eur_pol = 0.0, m_asia_pol = 0.0):
  demography = msprime.Demography()
  demography.add_population(name = "admix", initial_size = N_admix * math.exp(t1 * gr), growth_rate = gr)
  demography.add_population(name = "afr", initial_size = N_afr)
  demography.add_population(name = "eur", initial_size = N_eur)
  demography.add_population(name = "asia", initial_size = N_asia)
  demography.add_population(name = "pol", initial_size = N_pol)
  
  demography.add_admixture(time=t1, derived="admix", ancestral=["afr", "eur", "asia", "pol"], proportions = [r1, r2, r3, 1-r1-r2-r3])
  
  demography.add_population_split(time=t2, derived=["pol"], ancestral="asia")
  demography.add_population_parameters_change(time=t2, initial_size = N_aa, growth_rate=0, population="asia")
  
  demography.add_population_split(time=t3, derived=["asia"], ancestral="eur")
  demography.add_population_parameters_change(time=t3, initial_size = N_ooa, growth_rate=0, population="eur")
  
  demography.add_population_split(time=t4, derived=["eur"], ancestral="afr")
  demography.add_population_parameters_change(time=t4, initial_size = N_anc, growth_rate=0, population="afr")
  return demography


########## Latinos (the same demography as native hawaiians, but with different default parameters) ###########
def latinos_demo(t1 = 25, t2 = 253, t3 = 1018, t4 = 2094, r1 = 0.107, r2 = 0.442, r3 = 0.0, 
                 N_admix = 41579, N_afr = 4986, N_eur = 13341, N_asia = 10000, N_pol = 73170, N_aa = 3092, N_ooa = 2948, N_anc = 2846, gr = 0.132, 
                 m_afr_eur = 0.0, m_afr_asia = 0.0, m_afr_pol = 0.0, m_eur_asia = 0.0, m_eur_pol = 0.0, m_asia_pol = 0.0):
  Q = np.array([[-m_afr_eur-m_afr_asia-m_afr_pol, m_afr_eur, m_afr_asia, m_afr_pol], 
                [m_afr_eur, -m_afr_eur-m_eur_asia-m_eur_pol, m_eur_asia, m_eur_pol], 
                [m_afr_asia, m_eur_asia, -m_afr_asia-m_eur_asia-m_asia_pol, m_asia_pol],
                [m_afr_pol, m_eur_pol, m_asia_pol, -m_afr_pol-m_eur_pol-m_asia_pol]]) if max(m_afr_eur, m_afr_asia, m_afr_pol, m_eur_asia, m_eur_pol, m_asia_pol) > 0 else None
  
  demo = Demo()
  demo.add_phase(Phase(0, t1, [1/N_admix  * math.exp(-t1 * gr)], [gr], populations = ["admix"]))
  P_admixture = np.array([
      [r1, r2, r3, 1-r1-r2-r3]
  ])
  demo.add_phase(Phase(t1, t2, [1/N_afr, 1/N_eur, 1/N_asia, 1/N_pol], P = P_admixture, Q = Q, populations = ["afr", "eur", "asia", "pol"]), discretize = 100)
  P_pol_split = np.array([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
      [0, 0, 1]
  ])
  demo.add_phase(Phase(t2, t3, [1/N_afr, 1/N_eur, 1/N_aa], P = P_pol_split, populations = ["afr", "eur", "asia"]))
  P_asia_split = np.array([
      [1, 0],
      [0, 1],
      [0, 1]
  ])
  demo.add_phase(Phase(t3, t4, [1/N_afr, 1/N_ooa], P = P_asia_split, populations = ["afr", "eur"]))
  P_ooa_split = np.array([
      [1],
      [1]
  ])
  demo.add_phase(Phase(t4, math.inf, [1/N_anc], P = P_ooa_split))
  return demo

def latinos_demography(t1 = 25, t2 = 253, t3 = 1018, t4 = 2094, r1 = 0.107, r2 = 0.442, r3 = 0.0, 
                       N_admix = 41579, N_afr = 4986, N_eur = 13341, N_asia = 10000, N_pol = 73170, N_aa = 3092, N_ooa = 2948, N_anc = 2846, gr = 0.132, 
                       m_afr_eur = 0.0, m_afr_asia = 0.0, m_afr_pol = 0.0, m_eur_asia = 0.0, m_eur_pol = 0.0, m_asia_pol = 0.0):
  demography = msprime.Demography()
  demography.add_population(name = "admix", initial_size = N_admix * math.exp(t1 * gr), growth_rate = gr)
  demography.add_population(name = "afr", initial_size = N_afr)
  demography.add_population(name = "eur", initial_size = N_eur)
  demography.add_population(name = "asia", initial_size = N_asia)
  demography.add_population(name = "pol", initial_size = N_pol)
  
  demography.add_admixture(time=t1, derived="admix", ancestral=["afr", "eur", "asia", "pol"], proportions = [r1, r2, r3, 1-r1-r2-r3])
  
  demography.add_population_split(time=t2, derived=["pol"], ancestral="asia")
  demography.add_population_parameters_change(time=t2, initial_size = N_aa, growth_rate=0, population="asia")
  
  demography.add_population_split(time=t3, derived=["asia"], ancestral="eur")
  demography.add_population_parameters_change(time=t3, initial_size = N_ooa, growth_rate=0, population="eur")
  
  demography.add_population_split(time=t4, derived=["eur"], ancestral="afr")
  demography.add_population_parameters_change(time=t4, initial_size = N_anc, growth_rate=0, population="afr")
  return demography
