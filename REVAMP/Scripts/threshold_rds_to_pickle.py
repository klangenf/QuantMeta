#!/usr/bin/env python
"""
Convert R RDS RMSE threshold functions to Python pickle format.
"""

import pickle
from pathlib import Path

def create_func1_model():
    """Create func1 model (0-10 reads/bp) based on R model summary"""
    # From R summary: lm(formula = log(RMSE) + 0.25 ~ log(total_avg_depth))
    # Coefficients:
    # (Intercept)           0.89422
    # log(total_avg_depth)  0.68615

    return {
        'intercept': 0.89422,
        'slope': 0.68615,
        'type': 'log_linear'
    }

def create_func2_model():
    """Create func2 model (10-100 reads/bp) based on R model summary"""
    # From R summary: lm(formula = log(RMSE) + 0.25 ~ log(total_avg_depth))
    # Coefficients:
    # (Intercept)           0.88053
    # log(total_avg_depth)  0.77109

    return {
        'intercept': 0.88053,
        'slope': 0.77109,
        'type': 'log_linear'
    }

def create_func3_model():
    """Create func3 model (100-1000 reads/bp) based on R model summary"""
    # From R summary: lm(formula = log(RMSE) + 0.25 ~ log(total_avg_depth))
    # Coefficients:
    # (Intercept)           0.61270
    # log(total_avg_depth)  0.85463

    return {
        'intercept': 0.61270,
        'slope': 0.85463,
        'type': 'log_linear'
    }

def create_func4_model():
    """Create func4 model (>=1000 reads/bp) based on R model summary"""
    # Constant value: 2027.23

    return {
        'value': 2027.23,
        'type': 'constant'
    }

def main():
    """Convert all R RMSE threshold functions to Python pickle format"""
    output_dir = Path("/Users/kathrynlangenfeld/Documents/UMich_Postdoc/QuantMeta/REVAMP/Regressions/threshold_read_depth_variability")
    output_dir.mkdir(parents=True, exist_ok=True)

    models = {
        'RMSE_limit_func1': create_func1_model(),
        'RMSE_limit_func2': create_func2_model(),
        'RMSE_limit_func3': create_func3_model(),
        'RMSE_limit_func4': create_func4_model()
    }

    # Save each model as pickle
    for name, model in models.items():
        output_path = output_dir / f"{name}.pkl"
        with open(output_path, 'wb') as f:
            pickle.dump(model, f)
        print(f"Saved {name} to {output_path}")

    print("All RMSE threshold functions converted to pickle format!")

if __name__ == "__main__":
    main()