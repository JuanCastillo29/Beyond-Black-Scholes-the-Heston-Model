\# Beyond Black-Scholes: A Journey into Stochastic Volatility with the Heston Model and Real Market Calibration



In times of market stress, the assumptions behind Black-Scholes (constant volatility and log-normal price distribution) start to break down. This project is an exploration of how the Heston stochastic volatility model, combined with real option chain data and Monte Carlo simulation, offers a more flexible and accurate framework for option pricing, volatility surface fitting, and risk management.



This project investigates the performance of the Heston stochastic volatility model in pricing equity options. Using historical option chain data from `^SPX`, I calibrate the model to real market prices and simulate paths using Monte Carlo techniques. Results show improved fit to implied volatility surfaces and demonstrate the model’s advantages in capturing real market behaviors such as volatility clustering and skew. A delta-hedging simulation reveals the practical risks of using mis-specified models like Black-Scholes in a stochastic volatility environment.



---



\## Project Structure



* `Data/`: Raw and processed datasets.
* `Notebooks/`: Jupyter notebooks for data collection, modelling and backtesting.
* `Scripts/`: Reusable Python scripts
* `Outputs/`: Backtesting results, plots, performance metrics
* `requirements.txt`: Dependencies



---



\## Goals



---



\## Getting started



1\\. Clone the repository



```bash



git clone …



cd …



```



2\\. Create and activate a virtual environment (optional but recommended)



```bash



python3 -m venv venv



source venv/bin/activate  # On Windows: venv\\\\Scripts\\\\activate



```

3\\. Install dependencies listed in 'requirements.txt'



```bash 



pip install -r requirements.txt



```



4\\. Run notebooks in sequence for data collection and anaysis



---



References:





\\\* Project by Juan Castillo \*\\

