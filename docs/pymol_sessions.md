## PyMOL Sessions

A problem is that all operations in PyMOL need to be self-contained.
Multiple threads may be interective with PyMOL at the same time therefore it is important not to contaminate each other.

## Classic pymol with decorator
My initial solution was to implement a lock on PyMOL allowing only one singleton session. The old file is [here](original_transpiler.py).
To make sure the lock is released when an error happened, say parsing a `mYpRoTeIn.exe.pdb`, it was implemented as a decorator class that wrapped the bound methods which did a self contained PyMOL operation.

    class PyMolTranspilerDeco:
        """
        Decorator for the bound methods of PyMolTranspiler that use Pymol.
        The session is shared... so only one thread at the time ought to use PyMOL.
        If a session raises an error, it should be caught so everyhting is cleaned closed and the error raised for the logger.
        Conor has rightfully suggested that the lock should be handled by the scheduler. I.e. a request is made and the a job is added to a queue.
        Currently, each extra concurrent thread simply waits or dies if it waits too long.
        :var lock: the lock. A class attribute.
        :vartype lock: threading.Lock
        """
        lock = threading.Lock()
    
        def __init__(self, fun):
            functools.update_wrapper(self, fun)
            self.fun = fun
    
    
        def __get__(self, obj, type=None):
            return self.__class__(self.fun.__get__(obj, type))
    
        def __call__(self, *args, **kwargs):
            try:
                self.start_up()
                reply = self.fun(*args, **kwargs)
                self.close_up()
                return reply
            except Exception as err:
                if self.lock.locked(): #the code errored before the lock could be released.
                    self.close_up()
                    print('ERROR: ', err)
                raise err
    
        def clean_up(self):
            """
            Reset the Pymol instance without calling reintialise
            """
            pymol.cmd.remove('all')
            pymol.cmd.delete('all')
    
        def close_up(self):
            """
            Calls ``clean_up`` and releases the lock.
            """
            self.clean_up()
            if self.lock.locked():
                self.lock.release()
                PyMolTranspiler.current_task = f'[{datetime.utcnow()}] idle'
            else:
                warn('The lock was off already...')
    
        def start_up(self):
            """
            Starts the task in ``self.fun`` and takes the lock or waits.
            """
            if not self.lock.acquire(timeout=60): #one minute wait.
                self.clean_up() #something failed very very ungracefully.
                self.lock.acquire()
                PyMolTranspiler.current_task = f'[{datetime.utcnow()}] working.'
                warn('The thread waited for over a minute!')
            self.clean_up()

## pymol2.PyMOL()

Then I switched to `pymol2.PyMOL()`.
This is [not documented](https://github.com/schrodinger/pymol-open-source/tree/master/modules/pymol2) and has some glitches, which suggests it should be treated with cation.
Say, this will _kill_ the kernel:

    with pymol2.PyMOL() as p:
        p.stop()
        
It dies due to a segmentation fault (`exit(139)`/`SIGSEGV`(signal 11)). This is rather problematic as it brings down the whole server.
Catching this with `signal` would require the offending thread to be killed.

The context manager ought to close up everything gracefully even when there is a return or an exception in there (_cf._ [context manager documetation](https://docs.python.org/3/library/contextlib.html)).
My guess is that it 

So there are two things I needed to check: it does what it says on the tin (no leakage) and that it closes the sessions properly.

### Test 1

It does indeed do as it says on the tin.

    s = [pymol2.PyMOL() for i in range(10)]
    for i in range(10):
        s[i].start()
        assert s[i].cmd.select('*') == 0, 'ought to be empty'
        s[i].cmd.fetch('1ubq')
    for i in range(10):
        s[i].cmd.remove(f'resi {i+1}')
        assert s[i].cmd.select(f'resi {i+1}') == 0, 'ought to be removed just now'
    for i in range(1,10):
        assert s[i].cmd.select(f'resi {i} and name CA') == 1, 'ought to been there just now'
    for i in range(10):
        assert s[i].cmd.select(f'resi {i+1}') == 0, 'ought to be removed before'
        
No weirdness happened.

### Test 2

It does delete everything:

    import psutil
    for j in range(10):
        print(psutil.virtual_memory().used*1e-9)
        for i in range(10000):
            with pymol2.PyMOL() as p:
                p.cmd.load('1ubq.cif')
                
Which gives out:

    6.89176576
    6.947459072
    6.9888368640000005
    7.03531008
    7.07796992
    7.1336099840000005
    7.18815232
    7.22540544
    
Ubiquitin is about 100 kilobytes. So 100_000 of them is 10 GB.
There was a modest creep in used memory, which is probably just the garbage collector or the Jupyter notebook.

## Wizard

One issue is that the `wizard` does not work with `cmd2` (`cmd` attribute of `pymol2`), which I found out too late. Luckily, the solution is simple,
making a class for the contextmanager that uses Singleton and waits (cf. the two mutagenesis methods)

    class GlobalPyMOL(): #singleton but that waits for the other thread to release it.
        pymol = pymol2.SingletonPyMOL()
        pymol.start()
        pylock = Lock()
    
        def __init__(self):
            pass
    
        def __enter__(self):
            if not self.pylock.acquire(timeout=60):
                # something hung up.
                self.pymol.cmd.remove('*')
                self.pymol.cmd.delete('*')
                self.pylock.release() #pointless roundtrip to be safe.
                self.pylock.acquire()
                return self.pymol
    
            else:
                self.pymol.cmd.delete('*')
                return self.pymol
    
        def __exit__(self, exc_type, exc_val, exc_tb):
            self.pymol.cmd.delete('*')
            self.pylock.release()

## Huston, we have a problem...

In case of disfuctionality replace `import pymol2` with

    if __name__ == '__main__':
        pymol_argv = ['pymol', '-qc']
    else:
        import __main__
        __main__.pymol_argv = ['pymol', '-qc']
    import pymol
    pymol.finish_launching()
    pymol.cmd.set('fetch_path', os.getcwd() + '/michelanglo_app/temp')

and change all:

        with pymol2.PyMOL() as self.pymol:

to

        self.pymol = pymol
        if 1==1:
        
