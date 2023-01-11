import dataclasses


class DataClassUnpack:
    CACHE = {}

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
