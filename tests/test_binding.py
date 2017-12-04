# The myplugin module must be locatable by Python.
# If you configured CMake in the build directory ``/path/to/repo/build`` then,
# assuming you are in ``/path/to/repo``, run the tests with something like
#     PYTHONPATH=./build/src/pythonmodule pytest tests

def test_dependencies():
    import gmx
    assert gmx
    import gmx.core
    # holder = gmx.core.get_holder()
    # gmx.core.get_name();

def test_imports():
    import myplugin
    assert myplugin
    import gmx.core

def test_add_potential():
    import gmx
    import myplugin
    import pytest
    # gmx.data provides a sample minimal tpr file
    from gmx.data import tpr_filename
    system = gmx.System._from_file(tpr_filename)
    potential = myplugin.MyRestraint()
    generic_object = object()
    with pytest.raises(Exception) as exc_info:
        potential.bind(generic_object)
    assert str(exc_info).endswith("bind method requires a python capsule as input")

    system.add_potential(potential)
    with gmx.context.DefaultContext(system.workflow) as session:
        session.run()

def test_plugin_potential():
    import gmx
    import myplugin
    from gmx.data import tpr_filename
    system = gmx.System._from_file(tpr_filename)
    potential = myplugin.HarmonicRestraint()
    potential.set_params(1, 4, 2.0, 100.0)
    # potential.set_params(1, 4, 0, 0)

    system.add_potential(potential)
    with gmx.context.DefaultContext(system.workflow) as session:
        session.run()
