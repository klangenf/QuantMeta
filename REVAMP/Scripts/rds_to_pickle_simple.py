#!/usr/bin/env python
"""
Convert R RDS regression models to Python pickle format.
Simplified version without numpy dependency issues.
"""

import pickle
from pathlib import Path

def create_bin1_model():
    """Create bin1 model (0-10 reads/bp) based on R model summary"""
    # From R summary: lm(formula = log(avg_depth + 1) ~ poly(avg_GC, 2, raw = TRUE) +
    #                     poly(total_avg_depth, 2, raw = TRUE) + poly(E_rel, 2, raw = TRUE)
    coeffs = [
        4.960e+00,   # intercept
        -2.078e-01,  # avg_GC
        1.890e-01,   # avg_GC^2
        2.932e-01,   # total_avg_depth
        -1.217e-02,  # total_avg_depth^2
        -1.292e+01,  # E_rel
        8.596e+00    # E_rel^2
    ]

    return {
        'coefficients': coeffs,
        'response_type': 'log',
        'bin': 'bin1',
        'X_test': None,
        'y_original': None,
        'residuals': None,
        'rmse': None
    }

def create_bin2_model():
    """Create bin2 model (10-100 reads/bp) based on R model summary"""
    coeffs = [
        -0.676354,   # intercept
        0.230358,    # avg_GC
        -0.154324,   # avg_GC^2
        1.327062,    # log(total_avg_depth)
        -0.047577    # log(total_avg_depth)^2
    ]

    return {
        'coefficients': coeffs,
        'response_type': 'log',
        'bin': 'bin2',
        'X_test': None,
        'y_original': None,
        'residuals': None,
        'rmse': None
    }

def create_bin3_model():
    """Create bin3 model (100-1000 reads/bp) based on R model summary"""
    coeffs = [
        -0.0899075,  # intercept
        -0.1847100,  # avg_GC
        0.3711778,   # avg_GC^2
        1.0031344    # log(total_avg_depth)
    ]

    return {
        'coefficients': coeffs,
        'response_type': 'log',
        'bin': 'bin3',
        'X_test': None,
        'y_original': None,
        'residuals': None,
        'rmse': None
    }

def create_bin4_model():
    """Create bin4 model (>=1000 reads/bp) based on R model summary"""
    coeffs = [
        -22978.438,  # intercept
        -1745.171,   # avg_GC
        531.544,     # avg_GC^2
        3426.407     # log(total_avg_depth)
    ]

    return {
        'coefficients': coeffs,
        'response_type': 'linear',
        'bin': 'bin4',
        'X_test': None,
        'y_original': None,
        'residuals': None,
        'rmse': None
    }

def main():
    """Convert all R models to Python pickle format"""
    output_dir = Path("/Users/kathrynlangenfeld/Documents/UMich_Postdoc/QuantMeta/REVAMP/Regressions/read_depth_variability")
    output_dir.mkdir(parents=True, exist_ok=True)

    models = {
        'reg1': create_bin1_model(),
        'reg2': create_bin2_model(),
        'reg3': create_bin3_model(),
        'reg4': create_bin4_model()
    }

    # Save each model as pickle
    for name, model in models.items():
        output_path = output_dir / f"{name}.pkl"
        with open(output_path, 'wb') as f:
            pickle.dump(model, f)
        print(f"Saved {name} to {output_path}")

    print("All models converted to pickle format!")

if __name__ == "__main__":
    main()