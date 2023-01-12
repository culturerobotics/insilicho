import dataclasses
import typing

from insilicho import parameters


class DataClassUnpack:
    """A utility function to unpack dataclasses with defaults."""

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
