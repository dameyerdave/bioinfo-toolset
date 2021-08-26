from liftover import get_lifter


def liftover(_from, to, chromosome, position):
    converter = get_lifter(_from, to)
    target = converter[chromosome][position]
    if len(target) > 0:
        return target[0][0].replace('chr', ''), int(target[0][1])
    else:
        raise ValueError(
            f"No conversion results found for chr{chromosome}:{position} ({_from} -> {to})")
