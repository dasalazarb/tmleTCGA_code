import math
import random
from typing import List, Tuple

random.seed(42)


def sigmoid(x: float) -> float:
    # guard against overflow
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    else:
        z = math.exp(x)
        return z / (1.0 + z)


def dot(weights: List[float], features: List[float]) -> float:
    return sum(w * f for w, f in zip(weights, features))


def add_intercept(X: List[List[float]]) -> List[List[float]]:
    return [[1.0] + row for row in X]


def train_logistic(
    X: List[List[float]],
    y: List[int],
    lr: float = 0.05,
    epochs: int = 800,
) -> List[float]:
    weights = [0.0 for _ in range(len(X[0]))]
    for _ in range(epochs):
        gradient = [0.0 for _ in weights]
        for xi, yi in zip(X, y):
            pred = sigmoid(dot(weights, xi))
            error = pred - yi
            for j, xij in enumerate(xi):
                gradient[j] += error * xij
        n = float(len(X))
        weights = [w - lr * g / n for w, g in zip(weights, gradient)]
    return weights


def predict_logistic(weights: List[float], X: List[List[float]]) -> List[float]:
    return [sigmoid(dot(weights, xi)) for xi in X]


def generate_synthetic(n: int = 150) -> Tuple[List[dict], List[int], List[int]]:
    data = []
    Y, Delta, A = [], [], []
    for _ in range(n):
        # clinical features
        age = random.gauss(55, 8)
        grade = random.choice([0, 1])
        # omics features
        expr1 = random.gauss(0, 1)
        expr2 = random.gauss(0, 1)
        cnv1 = random.gauss(0, 1)

        # treatment assignment depends mostly on clinical vars
        p_treat = sigmoid(-1.0 + 0.04 * age + 0.8 * grade)
        a_val = 1 if random.random() < p_treat else 0

        # missingness indicator
        p_obs = sigmoid(1.5 - 0.5 * a_val + 0.2 * grade)
        delta_val = 1 if random.random() < p_obs else 0

        # outcome with treatment effect
        lin_y = -2.0 + 0.03 * age + 0.6 * grade + 0.9 * a_val + 0.4 * expr1
        p_y = sigmoid(lin_y)
        y_val = 1 if random.random() < p_y else 0

        A.append(a_val)
        Delta.append(delta_val)
        Y.append(y_val)
        data.append(
            {
                "age": age,
                "grade": grade,
                "expr1": expr1,
                "expr2": expr2,
                "cnv1": cnv1,
                "A": a_val,
            }
        )
    return data, Y, Delta


def to_matrix(data: List[dict], columns: List[str]) -> List[List[float]]:
    return [[row[col] for col in columns] for row in data]


def summarize(preds: List[float]) -> Tuple[float, float]:
    mean_val = sum(preds) / len(preds)
    # simple std dev
    var = sum((p - mean_val) ** 2 for p in preds) / len(preds)
    return mean_val, math.sqrt(var)


def main() -> None:
    data, Y, Delta = generate_synthetic()

    clinical_cols = ["age", "grade"]
    full_cols = clinical_cols + ["expr1", "expr2", "cnv1", "A"]

    clinical_X = add_intercept(to_matrix(data, clinical_cols))
    full_X = add_intercept(to_matrix(data, full_cols))

    # 1) Propensity score
    treat_weights = train_logistic(clinical_X, [row["A"] for row in data])
    g1W = predict_logistic(treat_weights, clinical_X)

    # 2) Missingness
    miss_weights = train_logistic(clinical_X, Delta)
    pDelta = predict_logistic(miss_weights, clinical_X)

    # 3) Outcome regression using observed data only
    observed_idx = [i for i, d in enumerate(Delta) if d == 1]
    obs_X = [full_X[i] for i in observed_idx]
    obs_Y = [Y[i] for i in observed_idx]
    outcome_weights = train_logistic(obs_X, obs_Y)

    # Potential outcomes
    full_X_a1 = [row[:-1] + [1.0] for row in full_X]
    full_X_a0 = [row[:-1] + [0.0] for row in full_X]

    q1W = predict_logistic(outcome_weights, full_X_a1)
    q0W = predict_logistic(outcome_weights, full_X_a0)

    ate = (sum(q1W) / len(q1W)) - (sum(q0W) / len(q0W))

    print("=== Early integration example (Python fallback) ===")
    print(f"Samples: {len(data)} | Observed outcomes: {len(observed_idx)}")
    ps_mean, ps_sd = summarize(g1W)
    print(f"Mean propensity score: {ps_mean:.3f} (sd={ps_sd:.3f})")
    pdelta_mean, _ = summarize(pDelta)
    print(f"Average observation probability: {pdelta_mean:.3f}")
    print(f"Estimated ATE (g-computation): {ate:.3f}")
    print("Preview of potential outcomes:")
    for i in range(3):
        print(
            f"  id={i}: q0W={q0W[i]:.3f}, q1W={q1W[i]:.3f}, g1W={g1W[i]:.3f}, Delta={Delta[i]}"
        )


if __name__ == "__main__":
    main()
