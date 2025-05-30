import unittest
import os
import sys

def run_all_tests():
    # Ensure the 'tests' directory is in the Python module search path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, current_dir)

    # Discover and run all tests in the 'tests' directory
    loader = unittest.TestLoader()
    suite = loader.discover(start_dir=current_dir, pattern="test_*.py")

    # Run the test suite
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Exit with appropriate status
    if result.wasSuccessful():
        print("\nAll tests passed successfully!")
    else:
        print("\nSome tests failed. Please check the output above.")

if __name__ == "__main__":
    run_all_tests()