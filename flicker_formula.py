import numpy as np
from scipy.optimize import fsolve

def calculate_time(D=2.2746e-8, l=0.08, x=0.05, target_value=0.8, terms=100):
    """Calculate time (t) for given parameters."""
    def C_x_t(D, t, l, x, terms):
        summation = sum((1 / (2 * n + 1)) * np.exp(-D * (2 * n + 1)**2 * np.pi**2 * t / l**2) *
                        np.sin((2 * n + 1) * np.pi * x / l) for n in range(terms))
        return 1 - (4 / np.pi) * summation

    def equation_to_solve(t, D, l, x, target_value, terms):
        return C_x_t(D, t, l, x, terms) - target_value

    initial_guess = 1000
    solution = fsolve(equation_to_solve, initial_guess, args=(D, l, x, target_value, terms))
    if solution[0] > 0:
        return solution[0]
    else:
        raise ValueError("Calculated time is not positive. Check input parameters.")

def calculate_diffusion(t=10800, l=0.08, x=0.05, target_value=0.8, terms=100):
    """Calculate diffusion coefficient (D) for given parameters."""
    def C_x_t(D, t, l, x, terms):
        summation = sum((1 / (2 * n + 1)) * np.exp(-D * (2 * n + 1)**2 * np.pi**2 * t / l**2) *
                        np.sin((2 * n + 1) * np.pi * x / l) for n in range(terms))
        return 1 - (4 / np.pi) * summation

    def equation_to_solve(D, t, l, x, target_value, terms):
        return C_x_t(D, t, l, x, terms) - target_value

    initial_guess = 1e-8
    solution = fsolve(equation_to_solve, initial_guess, args=(t, l, x, target_value, terms))
    if solution[0] > 0:
        return solution[0]
    else:
        raise ValueError("Calculated D is not positive. Check input parameters.")

def calculate_target(D=2.2746e-8, t=10800, l=0.08, x=0.05, terms=100):
    """Calculate target value (C(x, t)) for given parameters."""
    summation = sum((1 / (2 * n + 1)) * np.exp(-D * (2 * n + 1)**2 * np.pi**2 * t / l**2) *
                    np.sin((2 * n + 1) * np.pi * x / l) for n in range(terms))
    return 1 - (4 / np.pi) * summation

def main():
    while True:
        choice = input("""
What do you want to find?
    1. Time (t) - Enter 1
    2. Diffusion coefficient (D) - Enter 2
    3. Target value (C(x, t)) - Enter 3
    4. Exit - Enter any other key
Your choice: """)

        if choice == '1':
            D = float(input("Enter D (default 2.2746e-8 m^2/s): ") or 2.2746e-8)
            l = float(input("Enter l (default 0.08 m): ") or 0.08)
            x = float(input("Enter x (default 0.05 m): ") or 0.05)
            target_value = float(input("Enter target C(x, t) (default 0.8): ") or 0.8)
            t = calculate_time(D, l, x, target_value)
            print(f"Estimated time (t): {t:.2f} seconds")

        elif choice == '2':
            t = float(input("Enter t (default 10800 seconds): ") or 10800)
            l = float(input("Enter l (default 0.08 m): ") or 0.08)
            x = float(input("Enter x (default 0.05 m): ") or 0.05)
            target_value = float(input("Enter target C(x, t) (default 0.8): ") or 0.8)
            D = calculate_diffusion(t, l, x, target_value)
            print(f"Diffusion coefficient (D): {D:.6e} m^2/s")

        elif choice == '3':
            D = float(input("Enter D (default 2.2746e-8 m^2/s): ") or 2.2746e-8)
            t = float(input("Enter t (default 10800 seconds): ") or 10800)
            l = float(input("Enter l (default 0.08 m): ") or 0.08)
            x = float(input("Enter x (default 0.05 m): ") or 0.05)
            C_target = calculate_target(D, t, l, x)
            print(f"Target C(x, t): {C_target:.6f}")

        else:
            print("Exiting the program.")
            break

if __name__ == "__main__":
    main()