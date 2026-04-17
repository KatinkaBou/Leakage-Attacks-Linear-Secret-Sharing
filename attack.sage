# Author: Katharina Boudgoust
# Date: 2026-04-14
#
# Description:
#   SageMath implementation of distinguishing attacks on (additive and) Shamir secret sharing.
#
# Acknowledgement:
#   Parts of this code were developed with assistance from ChatGPT (OpenAI).

"""
==== HELPER FUNCTIONS ====
"""

def round_away_from_zero(x):
    """
    Round to nearest integer, with ties (.5) going away from 0.
    """
    if x >= 0:
        return Integer(floor(x + 1/2))
    else:
        return Integer(ceil(x - 1/2))
        
def centered_rep(x, p):
    """
    Return the centered representative of x in F_p, viewed in
    {-(p-1)//2, ..., (p-1)//2}.

    INPUT:
        - x : element of GF(p) or integer
        - p : odd prime

    OUTPUT:
        Integer representative of x in the centered interval.
    """
    if p % 2 == 0:
        raise ValueError("p must be odd")

    a = Integer(x) % p
    half = (p - 1) // 2

    if a <= half:
        return a
    else:
        return a - p   
        
def scaled_rounding(x, p, q):
    """
    Rounding map from F_p to Z:
      1. take centered representative in {-(p-1)//2, ..., (p-1)//2}
      2. scale by q/p
      3. round to nearest integer

    INPUT:
        - x : element of GF(p) or integer
        - p : odd prime
        - q : integer with 0 < q < p

    OUTPUT:
        Integer
    """
    if p % 2 == 0:
        raise ValueError("p must be odd")
    if not (0 < q < p):
        raise ValueError("q must satisfy 0 < q < p")

    x_centered = centered_rep(x, p)
    y = QQ(q * x_centered) / QQ(p) 
    return round_away_from_zero(y) 
    
def random_prime_of_bitlength(bits):
    """
    Return a random prime of exactly `bits` bits.

    INPUT:
        - bits : integer >= 2

    OUTPUT:
        - a prime integer p with exact bit size `bits`
    """
    if bits < 2:
        raise ValueError("bits must be at least 2")

    lower = 2**(bits - 1)
    upper = 2**bits - 1

    return random_prime(upper, lbound=lower) 
    
"""
==== SECRET SHARINGS ====
"""  
        
def additive_share(secret, n, p):
    """
    Generate an additive secret sharing of `secret` over F_p.

    INPUT:
        - secret : integer or element of GF(p)
        - n      : number of shares, integer >= 2
        - p      : prime modulus

    OUTPUT:
        - list of n shares in GF(p) whose sum is secret
    """
    if n < 2:
        raise ValueError("n must be at least 2")

    F = GF(p)
    s = F(secret)

    shares = [F.random_element() for _ in range(n - 1)]
    last_share = s - sum(shares, F(0))
    shares.append(last_share)
    return shares
    
def additive_reconstruct(shares, p):
    """
    Reconstruct the secret from additive shares over F_p.
    """
    F = GF(p)
    return sum((F(x) for x in shares), F(0))
    
def lagrange_coefficients_at_zero(xs, p):
    """
    Compute the Lagrange coefficients at 0 for the points xs over F_p.

    INPUT:
        - xs : list of distinct nonzero field elements (or integers mod p)
        - p  : odd prime

    OUTPUT:
        - list [lambda_1, ..., lambda_t] in GF(p) such that
              for any polynomial f of degree < len(xs),
              f(0) = sum_i lambda_i * f(xs[i])
    """
    F = GF(p)
    xs = [F(x) for x in xs]
    t = len(xs)

    if len(set(xs)) != t:
        raise ValueError("Interpolation points must be distinct")

    lambdas = []
    for i in range(t):
        num = F(1)
        den = F(1)
        xi = xs[i]
        for j in range(t):
            if j != i:
                xj = xs[j]
                num *= -xj
                den *= (xi - xj)
        lambdas.append(num / den)

    return lambdas

def shamir_share(secret, t, n, p, xs=None):
    """
    Generate a t-out-of-n Shamir secret sharing of `secret` over F_p.

    INPUT:
        - secret : element of F_p or integer mod p
        - t      : threshold, integer >= 1
        - n      : number of shares, integer >= t
        - p      : prime modulus
        - xs     : optional list of n distinct nonzero x-coordinates in F_p

    OUTPUT:
        - shares : list of pairs (x_i, y_i) in GF(p) x GF(p)
    """
    if t < 1:
        raise ValueError("t must be at least 1")
    if n < t:
        raise ValueError("n must be at least t")

    F = GF(p)
    secret = F(secret)

    if xs is None:
        xs = [F(i) for i in range(1, n + 1)]
    else:
        xs = [F(x) for x in xs]

    if len(xs) != n:
        raise ValueError("xs must have length n")
    if any(x == 0 for x in xs):
        raise ValueError("All x-coordinates must be nonzero")
    if len(set(xs)) != n:
        raise ValueError("x-coordinates must be distinct")

    coeffs = [secret] + [F.random_element() for _ in range(t - 1)]

    def f(x):
        return sum(coeffs[k] * (x ** k) for k in range(t))

    shares = [(x, f(x)) for x in xs]
    return shares  
    
def shamir_reconstruct_at_zero(shares, p):
    """
    Reconstruct the secret from a list of shares (x_i, y_i) using
    Lagrange interpolation at 0.
    """
    F = GF(p)
    xs = [F(x) for (x, _) in shares]
    ys = [F(y) for (_, y) in shares]
    lambdas = lagrange_coefficients_at_zero(xs, p)
    return sum(lambdas[i] * ys[i] for i in range(len(shares)))
    
def shamir_to_additive_shares(shares_subset, p):
    """
    Given a subset of Shamir shares of size t, transform them into additive
    shares z_i = lambda_i * y_i, where lambda_i are the Lagrange coefficients
    at 0 for the chosen x_i's.

    INPUT:
        - shares_subset : list of pairs (x_i, y_i)
        - p             : prime modulus

    OUTPUT:
        - additive_shares : list [z_1, ..., z_t] in GF(p)
    """
    F = GF(p)
    xs = [F(x) for (x, _) in shares_subset]
    ys = [F(y) for (_, y) in shares_subset]
    lambdas = lagrange_coefficients_at_zero(xs, p)
    additive_shares = [lambdas[i] * ys[i] for i in range(len(shares_subset))]
    return additive_shares
    
"""
==== DISTINGUISHING ATTACK WITH LEAKAGE ====
"""
    
def leakage_vector(shares, p, q):
    """
    Apply rounding to each share.
    Returns a list of integers.
    """
    return [scaled_rounding(share, p, q) for share in shares]

def challenge_messages(p):
    """
    Return the challenge messages (m0, m1) for the attack.

    For odd prime p:
      - m0 = 0
      - m1 = (p-1)/2 
    """
    if p % 2 == 0:
        raise ValueError("p must be odd")

    m0 = 0
    m1 = (p-1)//2
    
    return m0, m1
    
def challenger(m0, m1, t, n, p, q):
    """
    Security game challenger:
      - samples b in {0,1}
      - secret shares m_b
      - leaks rounded shares

    Returns:
      (b, shares, y)
    where y is the leakage vector.
    """
    b = ZZ.random_element(2)
    mb = m0 if b == 0 else m1
    shamir_shares = shamir_share(mb, t, n, p)
    chosen_shamir_shares =  [shamir_shares[i] for i in range(t)]
    shares = shamir_to_additive_shares(chosen_shamir_shares,p)
    #shares = additive_share(mb, n, p)
    y = leakage_vector(shares, p, q)
    return b, shares, y

def distinguisher(y, q, t):
    """
    Your distinguisher:
      - compute v = sum(y_i) mod q
      - take centered representative modulo q
      - guess 0 if |v| < t/2, else guess 1
    """
    v = Integer(sum(y)) % q
    #print("value v after modulo q: ",v)
    dist_to_0 = min(v, q - v)
    #print("distance to 0: ",dist_to_0)
    if dist_to_0 < QQ(t)/2:
        return 0
    else:
        return 1

def one_trial(m0, m1, t, n, p, q, verbose=False):
    """
    Run one challenge/distinguisher trial.
    Returns True iff the distinguisher guessed correctly.
    """
    b, shares, y = challenger(m0, m1, t, n, p, q)
    guess = distinguisher(y, q, t)
    
    if verbose:
        print(f"b      = {b}")
        print(f"shares = {shares}")
        print(f"y      = {y}")
        print(f"sum(y) = {sum(y)}")
        print(f"v mod q (centered) = {centered_rep(sum(y), q)}")
        print(f"guess  = {guess}")
        print(f"correct? {guess == b}")
    
    return guess == b

def experiment(bits, t, n, q, trials=1000, verbose=False):
    """
    Sample a random prime p of the given bit size, print it, and run the attack.

    INPUT:
        - bits    : bit length of p
        - t       : threshold number
        - n       : number of additive shares
        - q       : rounding modulus
        - trials  : number of trials
        - verbose : if True, print each trial too

    OUTPUT:
        dictionary containing p, p mod 4, messages, and success rate
    """
    p = random_prime_of_bitlength(bits)
    m0, m1 = challenge_messages(p)

    print("Selected prime p =", p)
    print("Bit length       =", p.nbits())
    print("m0               =", m0)
    print("m1               =", m1)
    print("q                =", q)
    print("n                =", n)
    print("t                =", t)
    print("trials           =", trials)

    wins = 0
    for _ in range(trials):
        if one_trial(m0, m1, t, n, p, q, verbose=False):
            wins += 1

    success_rate = QQ(wins) / QQ(trials)

    print("Success rate     =", success_rate)

    return {
        "p": p,
        "p_mod_4": p % 4,
        "m0": m0,
        "m1": m1,
        "q": q,
        "n": n,
        "trials": trials,
        "success_rate": success_rate
    }

print("")
print("=============")
print("Experiment 1")
print("=============")
print("")

bits = 10
q = 17
t = 4
n = 50

experiment(bits, t, n, q, trials=1000, verbose=False)

print("")
print("=============")
print("Experiment 2")
print("=============")
print("")

bits = 20
q = 17
t = 4
n = 50

experiment(bits, t, n, q, trials=1000, verbose=False)

print("")
print("=============")
print("Experiment 3")
print("=============")
print("")

bits = 40
q = 17
t = 4
n = 50

experiment(bits, t, n, q, trials=1000, verbose=False)

print("")
print("=============")
print("Experiment 4")
print("=============")
print("")

bits = 40
q = 17
t = 4
n = 100

experiment(bits, t, n, q, trials=1000, verbose=False)
