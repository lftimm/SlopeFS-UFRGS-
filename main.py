import SlopeFs as slf
import math
def main():
    """
    https://www.desmos.com/calculator/iibc0laiu1
    ^ Desmos pra verificar manualmente a geometria do problema.

    Como usar?
    SoilSpace() -> Define tanto o solo quanto o talude.
                Possui valores padrões que podem ser sobrescritos na criação da classe, ex: slope = SoilSpace(h=10)
                Os valores padrões estão na definição da classe lá em cima.
    SoilFS()    -> Envolve o código inteiro, com o modelo matemático, o cálculo do FS e sua plotagem.
                Você pode passar o talude como parâmetro para a inicialização da classe, caso não ele usará o SoilSpace() padrão.
                SoilFS.view() plota os resultados.
    """
    default = slf.SoilSpace()
    deg = math.degrees(math.atan(10/20))
    slope = slf.SoilSpace(h=10, gam=18.2, c=15, phi=20, alp=deg)

    fs = slf.SoilFs(slope, methods=['Fellenius', 'OSM'], minimize=True)
    #fs = SoilFs(default)

    print(fs)
    fs.view()


if __name__ == '__main__':
    main()
