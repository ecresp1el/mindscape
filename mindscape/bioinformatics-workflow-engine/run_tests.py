import unittest

def run_all_tests():
    # Create a test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add tests in the desired order
    suite.addTests(loader.loadTestsFromModule(__import__('tests.test_utils')))
    suite.addTests(loader.loadTestsFromModule(__import__('tests.test_base_workflow')))
    suite.addTests(loader.loadTestsFromModule(__import__('tests.test_workflow_manager')))

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