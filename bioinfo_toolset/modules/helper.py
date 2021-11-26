import json
import yaml


def filter_dict(d):
    try:
        return {k: filter_dict(d[k])
                for k in d
                if not k.startswith("_")}
    except Exception:
        return d


def as_str(s):
    return s if isinstance(s, str) else s.decode()


def dj(o):
    """print VR object as pretty formated json"""
    print(json.dumps(filter_dict(o.as_dict()), indent=2, sort_keys=True))


def dy(fns, o):
    """execute function f in fns on o, returning a yaml block representing the test"""
    r = {
        "in": o.as_dict(),
        "out": {f.__name__: as_str(f(o)) for f in fns}
    }
    print(yaml.dump(filter_dict(
        {o.type._value: {"-": r}})).replace("'-':", "-"))


