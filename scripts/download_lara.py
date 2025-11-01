import re
import urllib.request
from urllib.error import HTTPError, URLError

nuclides_str = '3H;6He;7Be;10Be;11C;13N;14C;15O;16N;18F;19O;22Na;23Ne;24Na;24Na-M;26Al;27Mg;28Mg;28Al;31Si;32P;33P;35S;36Cl;37S;37Ar;38Cl;40K;41Ar;41Ca;42K;43K;44Sc;44Ti;45Ca;45Ti;46Sc;47Ca;47Ca EQUI;47Sc;48Sc;48V;49Ca;50V;51Ti;51Cr;52V;52Mn;52Mn-M;52Fe EQUI;52Fe;54Mn;55Fe;55Co;56Mn;56Co;56Ni;57Co;57Ni;58Co-M;58Co;59Fe;59Ni;<option selected="">60Co;60Co-M;61Cu;63Ni;63Zn;64Cu;65Ni;65Zn;66Cu;66Ga;67Cu;67Ga;68Ga;68Ge EQUI;68Ge;69Zn;69Zn-M;72Ga;73Se;75Se;76As;76Br;77Ge;79Se;79Kr;80Br;81Kr-M;81Rb;81Rb-M;82Br;82Rb;82Sr;83Br;83Kr-M;83Rb;84Rb;85Kr;85Kr-M;85Sr;86Br;86Rb;87Br;87Kr;87Rb;88Kr;88Rb;88Y;88Zr;89Kr;89Rb;89Sr;89Y-M;89Zr;90Kr;90Rb-M;90Rb;90Sr;90Y;90Y-M;91Kr;91Rb;91Sr;91Y-M;91Y;92Sr;92Y;92Nb;92Nb-M;93Y;93Zr;93Nb-M;93Mo;93Mo-M;94Sr;94Y;94Nb;94Tc-M;95Zr;95Zr EQUI;95Nb;95Nb-M;96Nb;97Zr EQUI;97Zr;97Nb-M;97Nb;97Ru;99Mo;99Tc;99Tc-M;101Mo;101Tc;103Ru;103Rh-M;103Pd;104Tc;106Ru EQUI;106Ru;106Rh;108Ag-M;108Ag;109Pd;109Cd;110Ag;110Ag-M;111Ag;111In;113In-M;113Sn EQUI;113Sn;114In;114In-M EQUI;114In-M;115Cd-M;115Cd EQUI;115In-M;116In-M;117Sn-M;121Sn;121Sn-M;121Te;121Te-M;122Sb;123Te-M;123I;123Xe;124Sb;124I;125Sn;125Sb;125Sb EQUI;125Te-M;125I;125Xe;126Sb;126Sb-M;126I;127Sb;127Te;127Te-M;127Xe;128I;129Sn-M;129Te-M;129Te;129I;129Xe-M;131Te;131Te-M EQUI;131Te-M;131I;131Xe-M;131Cs;131Ba;132Sb;132Sb-M;132Te;132I;132I-M;132Xe-M;133Te;133Te-M EQUI;133Te-M;133I;133Xe;133Xe-M;133Ba;134Te;134I;134Cs;134Cs-M;135I EQUI;135I;135Xe-M;135Xe;135Cs;135Ba-M;136Cs;136Ba-M;137I;137Xe;137Cs;137Ba-M;138I;138Xe;138Cs;138La;139Xe;139Ba;139Ce;140Xe;140Cs;140Ba EQUI;140Ba;140La;141Ba;141La;141Ce;142Pr;143Cs;143Ce;143Pr;144Ce;144Ce EQUI;144Pr;144Pr-M;144Nd;146Pm;147Pr;147Nd;147Pm;147Sm;148Pm;148Pm-M;150Eu;151Pm;151Sm;152Eu-M;152Eu;153Sm;153Gd;154Eu;155Sm;155Eu;156Eu;159Gd;160Tb;161Tb;165Dy;166Ho;166Ho-M;169Er;169Yb-M;169Yb;170Tm;175Yb;175Hf;176Lu;176Lu-M;177Lu-M;177Lu;178Ta;178W;181Hf;181W;182Ta;182Ta-M;183Ta EQUI;183W;183W-M;184Re;184Re-M;185W;185W-M;185Os;186Re;187W;187Re;188W EQUI;188W;188Re;188Re-M;191Os;191Ir-M;191Pt;192Ir;193Os;193Pt;193Pt-M;194Ir;195Pt-M;195Au;196Au;197Pt;197Pt-M;197Au-M;197Hg-M;197Hg;198Au;199Au;200Tl;200Pb;201Tl;201Pb;202Tl;203Hg;203Pb;204Tl;205Pb;206Hg;206Tl;207Tl;207Bi;208Tl;209Tl;209Pb;209Bi;209Po;210Tl;210Pb;210Bi;210Po;211Pb;211Bi;211Po;211Po-M;211At;212Pb;212Bi;212Po;213Bi;213Bi EQUI;213Po;214Pb;214Bi;214Po;215Bi;215Po;215At;216Po;216At;216Rn;217At;217Rn;218Po;218At;218Rn;219At;219Rn;220Rn;220Fr;220Ra;221Fr;221Ra;222Rn;223Fr;223Ra;224Fr;224Ra;224Ac;225Ra;225Ac;226Ra EQUI;226Ra;227Fr;227Ra;227Ac;227Ac EQUI;227Th;228Ra;228Ac;228Th;228Pa;229Th;230Th;231Th;231Pa;232Th EQUI;232Th;232Pa;232U;233Th;233Pa;233U EQUI;233U;234Th;234Pa;234Pa-M;234U;235U EQUI;235U;235Np;236U;236Np-M;236Np;236Pu;237U;237Np;237Np EQUI;237Pu;238U;238U EQUI;238Np;238Pu;239U;239Np;239Pu;240U;240Np;240Pu;240Cm;241Pu;241Pu EQUI;241Am;242Pu;242Am-M;242Am;242Cm;243Pu;243Am;243Cm;244Pu;244Pu EQUI;244Am;244Am-M;244Cm;245Am;245Cm;246Cm;248Cm;252Cf'
nuclides_list = nuclides_str.split(";")


def check_url(url):
    try:
        urllib.request.urlopen(url)
        return True  # URL exists
    except HTTPError:
        return False  # Server returned an HTTP error
    except URLError:
        return False  # URL could not be reached


for this_nuclide in nuclides_list:
    match = re.match(r"([0-9]+)([a-z]+)", this_nuclide, re.I)
    if match:
        max_i = 0
        items = match.groups()
        mass_number = items[0]
        element = items[1]
        i = 0
        filename = f"{element}-{mass_number}_@{i:{0}>2}.lara.txt"
        url_str = "http://www.lnhb.fr/Laraweb/Results/" + filename
        if check_url(url_str):
            urllib.request.urlretrieve(url_str, filename)
        # for i in range(100):
        #    filename = f"{element}-{mass_number}_@{i:{0}>2}.lara.txt"
        #    url_str = "http://www.lnhb.fr/Laraweb/Results/" + filename
        #    print(f"Checking {url_str}")
        #    if check_url(url_str):
        #        print(f"{url_str} exists")
        #    else:
        #        print(f"{url_str} does not exist")
        #        break
        #    # urllib.request.urlretrieve(url_str, filename)
