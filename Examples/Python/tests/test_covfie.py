import acts
from acts import covfie_conversion as cc

def test_constant_field_conversion():
    v = acts.Vector3(1, 2, 3)
    af = acts.ConstantBField(v)
    cf = cc.covfieField(af)
    view = cc.newView(cf)
    for (x, y, z) in [(0,0,1), (1,1,1), (1,0,2)]:
        assert view.at(x, y, z) == [1, 2, 3]