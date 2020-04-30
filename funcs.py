

def tochka_peresecheniya(verhnaya_liniyax, verhnaya_liniyay, nizhnaya_liniyax2, nizhnaya_liniyay2):
    a_1 = verhnaya_liniyay[0] - verhnaya_liniyay[1]
    b_1 = verhnaya_liniyax[1] - verhnaya_liniyax[0]
    c_1 = verhnaya_liniyax[0] * verhnaya_liniyay[1] - verhnaya_liniyax[1] * verhnaya_liniyay[0]
    a_2 = nizhnaya_liniyay2[0] - nizhnaya_liniyay2[1]
    b_2 = nizhnaya_liniyax2[1] - nizhnaya_liniyax2[0]
    c_2 = nizhnaya_liniyax2[0] * nizhnaya_liniyay2[1] - nizhnaya_liniyax2[1] * nizhnaya_liniyay2[0]
    if a_1 * b_2 - a_2 * b_1 == 0:
        # print('Ошибка. Прямые не пересекаются')
        tochka_peresecheniyax = (verhnaya_liniyax[1] + nizhnaya_liniyax2[0]) / 2
        tochka_peresecheniyay = (verhnaya_liniyay[1] + nizhnaya_liniyay2[0]) / 2
        return [tochka_peresecheniyax, tochka_peresecheniyay]
    else:
        tochka_peresecheniyax = -(c_1 * b_2 - c_2 * b_1) / (a_1 * b_2 - a_2 * b_1)
        tochka_peresecheniyay = -(a_1 * c_2 - a_2 * c_1) / (a_1 * b_2 - a_2 * b_1)
        return [tochka_peresecheniyax, tochka_peresecheniyay]

