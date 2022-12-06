import msprime
from .glike import *

def threeway_admixture_demo(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e, m_ab, m_cd):
  demo = Demography()
  demo.add_phase(ContinuousPhase(0, t1, Q = np.zeros([1,1]), n = [1/N]))
  demo.add_phase(DiscretePhase(t1, P = np.array([[r1, 1-r1]])))
  demo.add_phase(ContinuousPhase(t1, t2, Q = np.array([[-m_ab, m_ab], [m_ab, -m_ab]]), n = [1/N_a, 1/N_b]))
  demo.add_phase(DiscretePhase(t2, P = np.array([[1, 0, 0], [0, r2, 1-r2]])))
  demo.add_phase(ContinuousPhase(t2, 100, Q = np.array([[0, 0, 0], [0, -m_cd, m_cd], [0, m_cd, -m_cd]]), n = [1/N_a, 1/N_c, 1/N_d]))
  demo.add_phase(ContinuousPhase(100, t3, Q = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]), n = [1/N_a, 1/N_c, 1/N_d]))
  demo.add_phase(DiscretePhase(t3, P = np.array([[1], [1], [1]])))
  demo.add_phase(ContinuousPhase(t3, 1e6, Q = np.array([[0]]), n = [1/N_e]))
  return demo

'''
def threeway_admixture_soft_demo(t1, t2, r1, r2, N, N_a, N_b, N_c, N_d, N_e, m_ab, m_cd):
  demo = Demography()
  demo.add_phase(ContinuousPhase(0, t1, Q = np.zeros([1,1]), n = [1/N]))
  demo.add_phase(DiscretePhase(t1, P = np.array([[r1, 1-r1]])))
  demo.add_phase(ContinuousPhase(t1, t2, Q = np.array([[-m_ab, m_ab], [m_ab, -m_ab]]), n = [1/N_a, 1/N_b]))
  demo.add_phase(DiscretePhase(t2, P = np.array([[1, 0, 0], [0, r2, 1-r2]])))
  demo.add_phase(ContinuousPhase(t2, 100, Q = np.array([[0, 0, 0], [0, -m_cd, m_cd], [0, m_cd, -m_cd]]), n = [1/N_a, 1/N_c, 1/N_d]))
  demo.add_phase(ContinuousPhase(100, 5e4, Q = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]), n = [1/N_a, 1/N_c, 1/N_d]))
  demo.add_phase(ContinuousPhase(5e4, 1e6, Q = np.array([[-1e-5, 1e-5, 1e-5], [1e-5, -1e-5, 1e-5], [1e-5, 1e-5, -1e-5]]), n = [1/N_a, 1/N_c, 1/N_d]))
  return demo
'''

def threeway_admixture_demography(t1, t2, r1, r2, N, N_a, N_b, N_c, N_d, N_e, m_ab, m_cd):
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
  
  demography.add_population_split(time=1e5, derived=["A", "C", "D"], ancestral="E")
  
  return demography

