from utils import CopelandGalaiCalc

model = CopelandGalaiCalc()

print("---Normal Distribution---")
model.bid_ask_normal_distribution(mu = 100, sigma = 10, pi= 0.3)

print("---Exponential Distribution---")
model.bid_ask_exponencial_distribution(lam = 0.5, pi= 0.3)



