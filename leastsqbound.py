"""Constrained multivariate least-squares optimization"""

import warnings

from numpy import array, take, eye, triu, transpose, dot, finfo
from numpy import empty_like, sqrt, cos, sin, arcsin, asarray
from numpy import atleast_1d, shape, issubdtype, dtype, inexact
from scipy.optimize import _minpack, leastsq


def _internal2external_grad(xi, bounds):
    """
    Calculate the internal (unconstrained) to external (constained)
    parameter gradiants.
    """
    grad = empty_like(xi)
    for i, (v, bound) in enumerate(zip(xi, bounds)):
        lower, upper = bound
        if lower is None and upper is None:  # No constraints
            grad[i] = 1.0
        elif upper is None:     # only lower bound
            grad[i] = v / sqrt(v * v + 1.)
        elif lower is None:     # only upper bound
            grad[i] = -v / sqrt(v * v + 1.)
        else:   # lower and upper bounds
            grad[i] = (upper - lower) * cos(v) / 2.
    return grad


def _internal2external_func(bounds):
    """
    Make a function which converts between internal (unconstrained) and
    external (constrained) parameters.
    """
    ls = [_internal2external_lambda(b) for b in bounds]

    def convert_i2e(xi):
        xe = empty_like(xi)
        xe[:] = [l(p) for l, p in zip(ls, xi)]
        return xe

    return convert_i2e


def _internal2external_lambda(bound):
    """
    Make a lambda function which converts a single internal (uncontrained)
    parameter to a external (constrained) parameter.
    """
    lower, upper = bound

    if lower is None and upper is None:  # no constraints
        return lambda x: x
    elif upper is None:     # only lower bound
        return lambda x: lower - 1. + sqrt(x * x + 1.)
    elif lower is None:     # only upper bound
        return lambda x: upper + 1. - sqrt(x * x + 1.)
    else:
        return lambda x: lower + ((upper - lower) / 2.) * (sin(x) + 1.)


def _external2internal_func(bounds):
    """
    Make a function which converts between external (constrained) and
    internal (unconstrained) parameters.
    """
    ls = [_external2internal_lambda(b) for b in bounds]

    def convert_e2i(xe):
        xi = empty_like(xe)
        xi[:] = [l(p) for l, p in zip(ls, xe)]
        return xi

    return convert_e2i


def _external2internal_lambda(bound):
    """
    Make a lambda function which converts an single external (constrained)
    parameter to a internal (unconstrained) parameter.
    """
    lower, upper = bound

    if lower is None and upper is None:  # no constraints
        return lambda x: x
    elif upper is None:     # only lower bound
        return lambda x: sqrt((x - lower + 1.) ** 2 - 1)
    elif lower is None:     # only upper bound
        return lambda x: sqrt((upper - x + 1.) ** 2 - 1)
    else:
        return lambda x: arcsin((2. * (x - lower) / (upper - lower)) - 1.)


def _check_func(checker, argname, thefunc, x0, args, numinputs,
                output_shape=None):
    res = atleast_1d(thefunc(*((x0[:numinputs],) + args)))
    if (output_shape is not None) and (shape(res) != output_shape):
        if (output_shape[0] != 1):
            if len(output_shape) > 1:
                if output_shape[1] == 1:
                    return shape(res)
            msg = "%s: there is a mismatch between the input and output " \
                  "shape of the '%s' argument" % (checker, argname)
            func_name = getattr(thefunc, '__name__', None)
            if func_name:
                msg += " '%s'." % func_name
            else:
                msg += "."
            raise TypeError(msg)
    if issubdtype(res.dtype, inexact):
        dt = res.dtype
    else:
        dt = dtype(float)
    return shape(res), dt


def leastsqbound(func, x0, args=(), bounds=None, Dfun=None, full_output=0,
                 col_deriv=0, ftol=1.49012e-8, xtol=1.49012e-8,
                 gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None):
    # use leastsq if no bounds are present
    if bounds is None:
        return leastsq(func, x0, args, Dfun, full_output, col_deriv,
                       ftol, xtol, gtol, maxfev, epsfcn, factor, diag)

    # create function which convert between internal and external parameters
    i2e = _internal2external_func(bounds)
    e2i = _external2internal_func(bounds)

    x0 = asarray(x0).flatten()
    i0 = e2i(x0)
    n = len(x0)
    if len(bounds) != n:
        raise ValueError('length of x0 != length of bounds')
    if not isinstance(args, tuple):
        args = (args,)
    shape, dtype = _check_func('leastsq', 'func', func, x0, args, n)
    m = shape[0]
    if n > m:
        raise TypeError('Improper input: N=%s must not exceed M=%s' % (n, m))
    if epsfcn is None:
        epsfcn = finfo(dtype).eps

    # define a wrapped func which accept internal parameters, converts them
    # to external parameters and calls func
    def wfunc(x, *args):
        return func(i2e(x), *args)

    if Dfun is None:
        if (maxfev == 0):
            maxfev = 200 * (n + 1)
        retval = _minpack._lmdif(wfunc, i0, args, full_output, ftol, xtol,
                                 gtol, maxfev, epsfcn, factor, diag)
    else:
        if col_deriv:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (n, m))
        else:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (m, n))
        if (maxfev == 0):
            maxfev = 100 * (n + 1)

        def wDfun(x, *args):  # wrapped Dfun
            return Dfun(i2e(x), *args)

        retval = _minpack._lmder(wfunc, wDfun, i0, args, full_output,
                                 col_deriv, ftol, xtol, gtol, maxfev,
                                 factor, diag)

    errors = {0: ["Improper input parameters.", TypeError],
              1: ["Both actual and predicted relative reductions "
                  "in the sum of squares\n  are at most %f" % ftol, None],
              2: ["The relative error between two consecutive "
                  "iterates is at most %f" % xtol, None],
              3: ["Both actual and predicted relative reductions in "
                  "the sum of squares\n  are at most %f and the "
                  "relative error between two consecutive "
                  "iterates is at \n  most %f" % (ftol, xtol), None],
              4: ["The cosine of the angle between func(x) and any "
                  "column of the\n  Jacobian is at most %f in "
                  "absolute value" % gtol, None],
              5: ["Number of calls to function has reached "
                  "maxfev = %d." % maxfev, ValueError],
              6: ["ftol=%f is too small, no further reduction "
                  "in the sum of squares\n  is possible.""" % ftol,
                  ValueError],
              7: ["xtol=%f is too small, no further improvement in "
                  "the approximate\n  solution is possible." % xtol,
                  ValueError],
              8: ["gtol=%f is too small, func(x) is orthogonal to the "
                  "columns of\n  the Jacobian to machine "
                  "precision." % gtol, ValueError],
              'unknown': ["Unknown error.", TypeError]}

    info = retval[-1]    # The FORTRAN return value

    if (info not in [1, 2, 3, 4] and not full_output):
        if info in [5, 6, 7, 8]:
            pass #warnings.warn(errors[info][0], RuntimeWarning)
        else:
            try:
                raise errors[info][1](errors[info][0])
            except KeyError:
                raise errors['unknown'][1](errors['unknown'][0])

    mesg = errors[info][0]
    x = i2e(retval[0])  # internal params to external params

    if full_output:
        # convert fjac from internal params to external
        grad = _internal2external_grad(retval[0], bounds)
        retval[1]['fjac'] = (retval[1]['fjac'].T / take(grad,
                             retval[1]['ipvt'] - 1)).T
        cov_x = None
        if info in [1, 2, 3, 4]:
            from numpy.dual import inv
            from numpy.linalg import LinAlgError
            perm = take(eye(n), retval[1]['ipvt'] - 1, 0)
            r = triu(transpose(retval[1]['fjac'])[:n, :])
            R = dot(r, perm)
            try:
                cov_x = inv(dot(transpose(R), R))
            except LinAlgError:
                pass
        return (x, cov_x) + retval[1:-1] + (mesg, info)
    else:
        return (x, info)
