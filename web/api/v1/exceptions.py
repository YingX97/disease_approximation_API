from flask import request
from flask_restful import abort



class FeatureStringFormatError(Exception):
    def __init__(self, msg, features):
        self.features = features
        super().__init__(self, msg)


def required_parameters(*required_args):
    """Decorator that aborts if mandatory parameters are missing."""
    def inner(wrapped):
        """Just an inner layer used by Python to get rid of decorator parameters."""
        def func(*args_inner, **kwargs_inner):
            """Decorated function."""
            for arg in required_args:
                if request.args.get(arg, None) is None:
                    abort(
                        400,
                        message=f"The \"{arg}\" parameter is required.",
                        error={
                            "type": "missing_parameter",
                            "missing_parameter": arg,
                        },
                    )
            return wrapped(*args_inner, **kwargs_inner)
        return func
    return inner


def model_exceptions(func):
    """Closure that deals with model exceptions at the API level."""
    def inner(*args_inner, **kwargs_inner):
        try:
            return func(*args_inner, **kwargs_inner)
        except Exception as exc:
            abort(
                500,
                message=exc
            )

    return inner