"""
Compare our X25519 implementation against the cryptography library.
"""

from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey, X25519PublicKey
from cryptography.hazmat.primitives import serialization

# Import our implementation
from ladder import *

# Curve25519 parameters
p = 2**255 - 19
r_bits = 256
BASE_POINT = 9

def test_against_cryptography_lib():
    print("=" * 70)
    print("Testing our X25519 against the cryptography library")
    print("=" * 70)
    
    # === Test 1: Public key generation ===
    print("\n[Test 1] Public key generation from private key")
    
    # Generate a keypair using the library
    lib_private_key = X25519PrivateKey.generate()
    lib_public_key = lib_private_key.public_key()
    
    # Get raw bytes
    private_bytes = lib_private_key.private_bytes_raw()
    public_bytes = lib_public_key.public_bytes_raw()
    
    # Convert to integers (little-endian)
    private_int = int.from_bytes(private_bytes, 'little')
    lib_public_int = int.from_bytes(public_bytes, 'little')
    
    print(f"  Private key (raw bytes): {private_bytes.hex()}")
    print(f"  Library public key:      {public_bytes.hex()}")
    
    # Compute public key with our implementation
    # Note: The library stores the unclamped private key, but clamps internally
    # We need to clamp it ourselves
    private_clamped = clamp_scalar(private_int)
    our_public_int = x25519_scalar_mult(private_clamped, BASE_POINT, p, r_bits)
    our_public_bytes = our_public_int.to_bytes(32, 'little')
    
    print(f"  Our public key:          {our_public_bytes.hex()}")
    
    if our_public_bytes == public_bytes:
        print("  ✓ Public keys match!")
    else:
        print("  ✗ Public keys DON'T match!")
        return False
    
    print_cycles()

    # === Test 2: Shared secret computation (ECDH) ===
    print("\n[Test 2] ECDH shared secret computation")
    
    # Generate two keypairs with the library
    alice_private = X25519PrivateKey.generate()
    alice_public = alice_private.public_key()
    
    bob_private = X25519PrivateKey.generate()
    bob_public = bob_private.public_key()
    
    # Library computes shared secret
    lib_shared_secret = alice_private.exchange(bob_public)
    
    print(f"  Alice private: {alice_private.private_bytes_raw().hex()}")
    print(f"  Alice public:  {alice_public.public_bytes_raw().hex()}")
    print(f"  Bob private:   {bob_private.private_bytes_raw().hex()}")
    print(f"  Bob public:    {bob_public.public_bytes_raw().hex()}")
    print(f"  Library shared secret: {lib_shared_secret.hex()}")
    
    # Our implementation computes shared secret
    alice_sk_int = int.from_bytes(alice_private.private_bytes_raw(), 'little')
    bob_pk_int = int.from_bytes(bob_public.public_bytes_raw(), 'little')
    
    alice_sk_clamped = clamp_scalar(alice_sk_int)
    bob_pk_masked = bob_pk_int & ((1 << 255) - 1)  # Clear top bit of public key
    
    our_shared_int = x25519_scalar_mult(alice_sk_clamped, bob_pk_masked, p, r_bits)
    our_shared_bytes = our_shared_int.to_bytes(32, 'little')
    
    print(f"  Our shared secret:     {our_shared_bytes.hex()}")
    
    if our_shared_bytes == lib_shared_secret:
        print("  ✓ Shared secrets match!")
    else:
        print("  ✗ Shared secrets DON'T match!")
        return False
    
    print_cycles()
    
    # === Test 3: Verify both sides compute same shared secret ===
    print("\n[Test 3] Verify ECDH symmetry (Alice and Bob get same secret)")
    
    # Bob computes shared secret using his private key and Alice's public key
    bob_sk_int = int.from_bytes(bob_private.private_bytes_raw(), 'little')
    alice_pk_int = int.from_bytes(alice_public.public_bytes_raw(), 'little')
    
    bob_sk_clamped = clamp_scalar(bob_sk_int)
    alice_pk_masked = alice_pk_int & ((1 << 255) - 1)
    
    bob_shared_int = x25519_scalar_mult(bob_sk_clamped, alice_pk_masked, p, r_bits)
    bob_shared_bytes = bob_shared_int.to_bytes(32, 'little')
    
    print(f"  Alice's computed secret: {our_shared_bytes.hex()}")
    print(f"  Bob's computed secret:   {bob_shared_bytes.hex()}")
    
    if our_shared_bytes == bob_shared_bytes:
        print("  ✓ Both sides compute the same secret!")
    else:
        print("  ✗ Secrets don't match!")
        return False
    
    # === Test 4: Multiple random tests ===
    print("\n[Test 4] Running 10 random ECDH exchanges...")
    
    for i in range(10):
        # Generate keypairs
        a_priv = X25519PrivateKey.generate()
        b_priv = X25519PrivateKey.generate()
        
        # Library shared secret
        lib_secret = a_priv.exchange(b_priv.public_key())
        
        # Our shared secret
        a_sk = clamp_scalar(int.from_bytes(a_priv.private_bytes_raw(), 'little'))
        b_pk = int.from_bytes(b_priv.public_key().public_bytes_raw(), 'little') & ((1 << 255) - 1)
        our_secret = x25519_scalar_mult(a_sk, b_pk, p, r_bits).to_bytes(32, 'little')
        
        if lib_secret != our_secret:
            print(f"  ✗ Test {i+1} failed!")
            return False
        print(f"  ✓ Test {i+1} passed")
    
    print("\n" + "=" * 70)
    print("All tests passed! Our X25519 implementation is correct.")
    print("=" * 70)
    return True

if __name__ == "__main__":
    test_against_cryptography_lib()
