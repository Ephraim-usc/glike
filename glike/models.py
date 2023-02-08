import msprime
from .glike import *


########## Twoway Admixture ###########
def twoway_admixture_demo(t1, t2, r, N, N_a, N_b, N_c):
  demo = Demo()
  demo.add_phase(Phase(0, [1/N]))
  demo.add_phase(Phase(t1, [1/N_a, 1/N_b], P = np.array([[r, 1-r]])))
  demo.add_phase(Phase(t2, [1/N_c], P = np.array([[1], [1]])))
  return demo

def twoway_admixture_demography(t1, t2, r, N, N_a, N_b, N_c, m):
  demography = msprime.Demography()
  demography.add_population(name = "O", initial_size = N)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.add_population(name = "C", initial_size = N_c)
  
  demography.add_admixture(time=t1, derived="O", ancestral=["A", "B"], proportions = [r, 1-r])
  demography.add_symmetric_migration_rate_change(time=t1, populations = ["A", "B"], rate = m)
  demography.add_symmetric_migration_rate_change(time=t2, populations = ["A", "B"], rate = 0)
  
  demography.add_population_split(time=t2, derived=["A", "B"], ancestral="C")
  
  return demography


########## Threeway Admixture ###########
def threeway_admixture_demo(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e):
  demo = Demo()
  demo.add_phase(Phase(0, [1/N]))
  demo.add_phase(Phase(t1, [1/N_a, 1/N_b], P = np.array([[r1, 1-r1]])))
  demo.add_phase(Phase(t2, [1/N_a, 1/N_c, 1/N_d], P = np.array([[1, 0, 0], [0, r2, 1-r2]])))
  demo.add_phase(Phase(t3, [1/N_e], P = np.array([[1], [1], [1]])))
  return demo

def threeway_admixture_demography(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e, m_ab, m_cd):
  demography = msprime.Demography()
  demography.add_population(name = "O", initial_size = N)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.add_population(name = "C", initial_size = N_c)
  demography.add_population(name = "D", initial_size = N_d)
  demography.add_population(name = "E", initial_size = N_e)
  
  demography.add_admixture(time=t1, derived="O", ancestral=["A", "B"], proportions = [r1, 1-r1])
  demography.add_symmetric_migration_rate_change(time=t1, populations = ["A", "B"], rate = m_ab)
  demography.add_symmetric_migration_rate_change(time=t2, populations = ["A", "B"], rate = 0)
  
  demography.add_admixture(time=t2, derived="B", ancestral=["C", "D"], proportions = [r2, 1-r2])
  demography.add_symmetric_migration_rate_change(time=t2, populations = ["C", "D"], rate = m_cd)
  demography.add_symmetric_migration_rate_change(time=100, populations = ["C", "D"], rate = 0)
  
  demography.add_population_split(time=t3, derived=["A", "C", "D"], ancestral="E")
  
  return demography


########## American Admixture ###########
def american_admixture_demo(t1, t2, t3, t4, r1, r2, N_afr, N_eur, N_asia, N_admix, N_ooa, N_anc, gr_eur, gr_asia, gr_admix):
  demo = Demo()
  demo.add_phase(Phase(0, [1/N_afr, (1/N_eur, gr_eur), (1/N_asia, gr_asia), (1/N_admix, gr_admix)],
                      populations = ["afr", "eur", "asia", "admix"]))
  P_admixture = np.array([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
      [r1, r2, 1-r1-r2]
  ])
  demo.add_phase(Phase(t1, [1/N_afr, (1/N_eur*math.exp(gr_eur*t1), gr_eur), (1/N_asia*math.exp(gr_asia*t1), gr_asia)], P = P_admixture,
                      populations = ["afr", "eur", "asia"]))
  P_asia_split = np.array([
      [1, 0],
      [0, 1],
      [0, 1]
  ])
  demo.add_phase(Phase(t2, [1/N_afr, 1/N_ooa], P = P_asia_split,
                      populations = ["afr", "eur"]))
  P_ooa_split = np.array([
      [1],
      [1]
  ])
  demo.add_phase(Phase(t3, [1/N_afr], P = P_ooa_split,
                       populations = ["afr"]))
  demo.add_phase(Phase(t4, [1/N_anc],
                       populations = ["afr"]))
  return demo

def american_admixture_demography(t1, t2, t3, t4, r1, r2, N_afr, N_eur, N_asia, N_admix, N_ooa, N_anc, gr_eur, gr_asia, gr_admix):
  demography = msprime.Demography()
  demography.add_population(name = "afr", initial_size = N_afr)
  demography.add_population(name = "eur", initial_size = N_eur, growth_rate = gr_eur)
  demography.add_population(name = "asia", initial_size = N_asia, growth_rate = gr_asia)
  demography.add_population(name = "admix", initial_size = N_admix, growth_rate = gr_admix)
  
  demography.add_admixture(time=t1, derived="admix", ancestral=["afr", "eur", "asia"], proportions = [r1, r2, 1-r1-r2])
  
  demography.add_mass_migration(time=t2, source="asia", dest="eur", proportion=1)
  demography.add_population_parameters_change(time=t2, initial_size = N_ooa, growth_rate=0, population="eur")
  
  demography.add_mass_migration(time=t3, source="eur", dest="afr", proportion=1)
  
  demography.add_population_parameters_change(time=t4, initial_size = N_anc, growth_rate=0, population="afr")
  return demography


########## Ancient Europe ###########
def ancient_europe_demo(t1, t2, t3, t4, t5, t6, r1, r2, r3, N_ana, N_neo, N_whg, N_bronze, N_yam, N_ehg, N_chg, N_ne, N_wa, N_ooa, gr):
  demo = Demo()
  demo.add_phase(Phase(0, [1/N_ana, 1/N_neo, 1/N_whg, (1/(N_bronze*math.exp(gr*t1)) , gr), 1/N_yam, 1/N_ehg, 1/N_chg],
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
  demo.add_phase(Phase(t1, [1/N_ana, 1/N_neo, 1/N_whg, 1/N_yam, 1/N_ehg, 1/N_chg], P = P_bronze_admixture,
                      populations = ["ana", "neo", "whg", "yam", "ehg", "chg"]))
  P_yam_admixture = np.array([
      [1, 0, 0, 0, 0],
      [0, 1, 0, 0, 0],
      [0, 0, 1, 0, 0],
      [0, 0, 0, r2, 1-r2],
      [0, 0, 0, 1, 0],
      [0, 0, 0, 0, 1]
  ])
  demo.add_phase(Phase(t2, [1/N_ana, 1/N_neo, 1/N_whg, 1/N_ehg, 1/N_chg], P = P_yam_admixture,
                      populations = ["ana", "neo", "whg", "ehg", "chg"]))
  P_neo_admixture = np.array([
      [1, 0, 0, 0],
      [r3, 1-r3, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1]
  ])
  demo.add_phase(Phase(t3, [1/N_ana, 1/N_whg, 1/N_ehg, 1/N_chg], P = P_neo_admixture,
                      populations = ["ana", "whg", "ehg", "chg"]))
  P_ne_split = np.array([
      [1, 0, 0],
      [0, 1, 0],
      [0, 1, 0],
      [0, 0, 1]
  ])
  demo.add_phase(Phase(t4, [1/N_ana, 1/N_ne, 1/N_chg], P = P_ne_split,
                      populations = ["ana", "ne", "chg"]))
  P_wa_split = np.array([
      [1, 0],
      [0, 1],
      [1, 0]
  ])
  demo.add_phase(Phase(t5, [1/N_wa, 1/N_ne], P = P_wa_split,
                      populations = ["wa", "ne"]))
  P_ooa_split = np.array([
      [1],
      [1]
  ])
  demo.add_phase(Phase(t6, [1/N_ooa], P = P_ooa_split,
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

