import numpy as np
from scipy.optimize import fsolve
import tkinter as tk
from tkinter import messagebox

def calculate_time(D, l, x, target_value, terms=100):
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

def calculate_diffusion(t, l, x, target_value, terms=100):
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

def calculate_target(D, t, l, x, terms=100):
    """Calculate target value (C(x, t)) for given parameters."""
    summation = sum((1 / (2 * n + 1)) * np.exp(-D * (2 * n + 1)**2 * np.pi**2 * t / l**2) *
                    np.sin((2 * n + 1) * np.pi * x / l) for n in range(terms))
    return 1 - (4 / np.pi) * summation

def on_calculate():
    try:
        choice = choice_var.get()
        if choice == 1:
            D = float(D_entry.get())
            l = float(l_entry.get())
            x = float(x_entry.get())
            target_value = float(target_entry.get())
            t = calculate_time(D, l, x, target_value)
            messagebox.showinfo("Result", f"Estimated time (t): {t:.2f} seconds")
        elif choice == 2:
            t = float(t_entry.get())
            l = float(l_entry.get())
            x = float(x_entry.get())
            target_value = float(target_entry.get())
            D = calculate_diffusion(t, l, x, target_value)
            messagebox.showinfo("Result", f"Diffusion coefficient (D): {D:.6e} m^2/s")
        elif choice == 3:
            D = float(D_entry.get())
            t = float(t_entry.get())
            l = float(l_entry.get())
            x = float(x_entry.get())
            C_target = calculate_target(D, t, l, x)
            messagebox.showinfo("Result", f"Target C(x, t): {C_target:.6f}")
    except Exception as e:
        messagebox.showerror("Error", str(e))

# Tkinter UI
root = tk.Tk()
root.title("Diffusion Calculator")

choice_var = tk.IntVar()
tk.Radiobutton(root, text="Calculate Time (t)", variable=choice_var, value=1).grid(row=0, column=0, sticky="w")
tk.Radiobutton(root, text="Calculate Diffusion (D)", variable=choice_var, value=2).grid(row=1, column=0, sticky="w")
tk.Radiobutton(root, text="Calculate Target (C(x, t))", variable=choice_var, value=3).grid(row=2, column=0, sticky="w")

tk.Label(root, text="D (m^2/s):").grid(row=3, column=0)
D_entry = tk.Entry(root)
D_entry.grid(row=3, column=1)
D_entry.insert(0, "2.2746e-8")

tk.Label(root, text="t (seconds):").grid(row=4, column=0)
t_entry = tk.Entry(root)
t_entry.grid(row=4, column=1)
t_entry.insert(0, "10800")

tk.Label(root, text="l (meters):").grid(row=5, column=0)
l_entry = tk.Entry(root)
l_entry.grid(row=5, column=1)
l_entry.insert(0, "0.08")

tk.Label(root, text="x (meters):").grid(row=6, column=0)
x_entry = tk.Entry(root)
x_entry.grid(row=6, column=1)
x_entry.insert(0, "0.05")

tk.Label(root, text="Target C(x, t):").grid(row=7, column=0)
target_entry = tk.Entry(root)
target_entry.grid(row=7, column=1)
target_entry.insert(0, "0.8")

tk.Button(root, text="Calculate", command=on_calculate).grid(row=8, column=0, columnspan=2)

root.mainloop()