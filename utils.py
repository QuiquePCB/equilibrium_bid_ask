import numpy as np
from scipy.stats import norm
from scipy.optimize import fsolve

class CopelandGalaiCalc: 

    def __init__(self):
        # Empty constructor
        pass
        
        
    def bid_ask_normal_distribution(self, mu:float=100, sigma:float=10, pi:float=0.3) -> tuple:
        """
        Calculates the bid and ask prices assuming the asset value follows
        a normal distribution.

        Parameters:
        - mu: Expected value of the asset
        - sigma: Standard deviation of the asset's value
        - pi: Probability that the trader is informed

        Returns:
        - A tuple with the ask and bid prices
        """

        def expected_profit(S):
            a = mu + S / 2 # Ask price
            b = mu - S / 2 # Bid price

            # Standarized normal variables
            z_a = (a - mu) / sigma
            z_b = (mu - b) / sigma

            # Conditional expected values for informed traders
            E_V_gt_a = mu + sigma * norm.pdf(z_a) / (1 - norm.cdf(z_a))
            E_V_lt_b = mu - sigma * norm.pdf(z_b) / norm.cdf(z_b)

            # Expected profits for each type of trader
            profit_ask = (1 - pi)*(a - mu) + pi*(a - E_V_gt_a)
            profit_bid = (1 - pi)*(mu - b) + pi*(E_V_lt_b - b)

            # Total expected profit
            expected_profit = 0.5 * (profit_ask + profit_bid)
            return expected_profit

        # Use fsolve to find the spread S that makes expected profit = 0
        initial_guess = 1.0
        spread_equilibrium = fsolve(expected_profit, initial_guess)[0]

        # Calculate the final ask and bid prices
        ask = mu + spread_equilibrium / 2
        bid = mu - spread_equilibrium / 2

        print(f"Equilibrium Spread: {spread_equilibrium:.4f}")
        print(f"Ask Price: {ask:.4f}")
        print(f"Bid Price: {bid:.4f}")

        return (ask, bid)
        
    def bid_ask_exponencial_distribution(self, lam: float=0.5, pi:float=0.3) -> tuple:
        """
        Calculates the bid and ask prices assuming the asset value follows
        an exponential distribution.

        Parameters:
        - lam: Lambda parameter of the exponential distribution
        - pi: Probability that the trader is informed

        Returns:
        - A tuple with the ask and bid prices
        """

        Ev = 1 / lam  # Expected value of the exponential distribution

        def expected_profit(S):

            a = Ev + (S / 2) # Ask price
            b = Ev - (S / 2) # Bid price

            # Conditional expected values for informed traders
            E_V_gt_a = a + 1/lam
            E_V_lt_b = (1/lam)-(b* np.exp(-b * lam)/1- np.exp(-lam*b))

            # Expected profits for each type of trader
            profit_ask = (1 - pi) * (a - Ev) + pi * (a - E_V_gt_a)
            profit_bid = (1 - pi) * (Ev - b) + pi * ( E_V_lt_b- b)

            return 0.5 * (profit_ask + profit_bid)

        # Solve for the spread S that makes expected profit = 0
        initial_guess = 1.0
        spread_equilibrium = fsolve(expected_profit, initial_guess)[0]

        # Final ask and bid prices
        ask = Ev + spread_equilibrium / 2
        bid = Ev - spread_equilibrium / 2

        print(f"Equilibrium Spread: {spread_equilibrium:.4f}")
        print(f"Ask Price: {ask:.4f}")
        print(f"Bid Price: {bid:.4f}")

        return (ask, bid)