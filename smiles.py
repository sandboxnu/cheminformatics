smiles = [{'smile': "Fake Smile"}, {}]

##or
class Smile(dict):
    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)

##
smile = dict([('smile', 'jeje'), ('murcko-smile', 'j')])