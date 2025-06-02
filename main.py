from utils import CopelandGalaiCalc

model = CopelandGalaiCalc()

mu = 480
sigma = 20
pi = 0.4

model.bid_ask_normal_distribution(mu, sigma, pi)

