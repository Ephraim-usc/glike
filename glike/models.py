import msprime
from .glike import *

########## American Admixture ###########
def american_admixture_demo(t1, t2, t3, t4, r1, r2, N_afr, N_eur, N_asia, N_admix, N_ooa, N_anc, gr_eur, gr_asia, gr_admix):
  demo = Demo()
  demo.add_phase(Phase(0, [1/N_afr, (1/N_eur, gr_eur), (1/N_asia, gr_asia), (1/N_admix, gr_admix)]))
  P_admixture = np.array([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
      [r1, r2, 1-r1-r2]
  ])
  demo.add_phase(Phase(t1, [1/N_afr, (1/N_eur*math.exp(gr_eur*t1), gr_eur), (1/N_asia*math.exp(gr_asia*t1), gr_asia)], P = P_admixture))
  P_asia_split = np.array([
      [1, 0],
      [0, 1],
      [0, 1]
  ])
  demo.add_phase(Phase(t2, [1/N_afr, 1/N_ooa], P = P_asia_split))
  P_ooa_split = np.array([
      [1],
      [1]
  ])
  demo.add_phase(Phase(t3, [1/N_afr], P = P_ooa_split))
  demo.add_phase(Phase(t4, [1/N_anc]))
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
  demo.add_phase(Phase(0, [1/N_ana, 1/N_neo, 1/N_whg, (1/N_bronze, gr), 1/N_yam, 1/N_ehg, 1/N_chg]))
  P_bronze_admixture = np.array([
      [1, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [0, 0, 1, 0, 0, 0],
      [0, r1, 0, 1-r1, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1]
  ])
  demo.add_phase(Phase(t1, [1/N_ana, 1/N_neo, 1/N_whg, 1/N_yam, 1/N_ehg, 1/N_chg], P = P_bronze_admixture))
  P_yam_admixture = np.array([
      [1, 0, 0, 0, 0],
      [0, 1, 0, 0, 0],
      [0, 0, 1, 0, 0],
      [0, 0, 0, r2, 1-r2],
      [0, 0, 0, 1, 0],
      [0, 0, 0, 0, 1]
  ])
  demo.add_phase(Phase(t2, [1/N_ana, 1/N_neo, 1/N_whg, 1/N_ehg, 1/N_chg], P = P_yam_admixture))
  P_neo_admixture = np.array([
      [1, 0, 0, 0],
      [r3, 1-r3, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1]
  ])
  demo.add_phase(Phase(t3, [1/N_ana, 1/N_whg, 1/N_ehg, 1/N_chg], P = P_neo_admixture))
  P_ne_split = np.array([
      [1, 0, 0],
      [0, 1, 0],
      [0, 1, 0],
      [0, 0, 1]
  ])
  demo.add_phase(Phase(t4, [1/N_ana, 1/N_ne, 1/N_chg], P = P_ne_split))
  P_wa_split = np.array([
      [1, 0],
      [0, 1],
      [1, 0]
  ])
  demo.add_phase(Phase(t5, [1/N_wa, 1/N_ne], P = P_wa_split))
  P_ooa_split = np.array([
      [1],
      [1]
  ])
  demo.add_phase(Phase(t6, [1/N_ooa], P = P_ooa_split))
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

