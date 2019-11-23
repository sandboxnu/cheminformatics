
def color_hex_to_array(hex):
    hex = hex.lstrip('#')

    return tuple(int(hex[i:i+2], 16) for i in (0, 2, 4))
