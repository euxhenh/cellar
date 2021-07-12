class InvalidArgument(Exception):
    pass


class InappropriateArgument(Exception):
    pass


class MethodNotImplementedError(NotImplementedError):
    pass


class IncorrectFileFormat(RuntimeError):
    pass


class InternalError(RuntimeError):
    pass


class UserError(RuntimeError):
    pass


class UIError(RuntimeError):
    pass
