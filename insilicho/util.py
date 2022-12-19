import dataclasses
import typing

from insilicho import parameters


class DataClassUnpack:
    CACHE: typing.Dict[
        typing.Union[parameters.InputParameters, parameters.InitialConditions],
        typing.Set[str],
    ] = {}

    @classmethod
    def instantiate(cls, cls_inst, arg_dict=None):
        if cls_inst not in cls.CACHE:
            cls.CACHE[cls_inst] = {
                f.name for f in dataclasses.fields(cls_inst) if f.init
            }

        if not arg_dict:
            return None

        return cls_inst(
            **{k: v for k, v in arg_dict.items() if k in cls.CACHE[cls_inst]}
        )


def sigmoid(x: typing.Any, sharpness: float = 1.0, thresh: float = 0.0):
    """

    Args:
        x: value
        sharpness : sharpness of transition. Defaults to 1.
        thresh : critical value for switching where function value is 0.5. Defaults to 0.
    """
    return 1 / (1 + np.exp(-sharpness * (x - thresh)))
