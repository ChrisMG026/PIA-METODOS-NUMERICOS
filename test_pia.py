import math
import unittest

import PIA


class PruebasMetodosNumericos(unittest.TestCase):
    def setUp(self):
        identificaciones = ("1234567", "2345678", "3456789", "4567890")
        self.matriculas = {
            f"M{indice}": [int(digito) for digito in texto]
            for indice, texto in enumerate(sorted(identificaciones), start=1)
        }

    def test_pendiente_global_usa_los_datos_del_documento(self):
        datos = PIA.calcular_documento_d()
        self.assertAlmostEqual(datos["k"], 1.71, places=12)
        self.assertEqual(datos["r_vals"], (1, 3, 5, 7))
        self.assertEqual(datos["y_vals"], (2.0, 5.3, 8.9, 12.2))

    def test_biseccion_selecciona_raiz_interior_del_documento_d(self):
        datos = PIA.calcular_documento_d_completo()
        self.assertAlmostEqual(datos["n"], 2.0, places=8)
        self.assertGreater(datos["iteraciones_n"], 0)
        self.assertAlmostEqual(datos["historial_n"][-1]["fx"], 0.0, delta=1e-8)

    def test_aislamiento_cubico_encuentra_tres_raices(self):
        raices = PIA.raices_reales_cubica(1, -6, 11, -6)
        self.assertEqual(len(raices), 3)
        for observada, esperada in zip(raices, (1, 2, 3)):
            self.assertAlmostEqual(observada, esperada, places=7)

    def test_f_y_g_son_finitos_y_coherentes_con_el_modelo(self):
        modelo = PIA.calcular_modelo(self.matriculas)
        variables = modelo["variables_a"]
        thermal = modelo["documento_f"]
        saturation = modelo["documento_c"]
        global_index = modelo["documento_g"]

        self.assertEqual(thermal["theta0"], 25.0)
        self.assertAlmostEqual(
            thermal["xi"],
            variables["a"]["valor"] + variables["d"]["valor"] + variables["g"]["valor"],
        )
        self.assertGreater(thermal["escala_fuente"], 1.0)
        self.assertAlmostEqual(
            thermal["fuente_maxima"] / thermal["escala_fuente"],
            thermal["fraccion_capacidad"] * thermal["capacidad_maxima"],
        )
        self.assertTrue(math.isfinite(thermal["theta_final"]))
        residual = (
            saturation["a_s"] * global_index["U_star"]**3
            - saturation["b_s"] * global_index["U_star"]**2
            + saturation["c_s"] * global_index["U_star"]
            - saturation["d_s"]
        )
        self.assertAlmostEqual(residual, 0.0, delta=1e-6)
        self.assertTrue(math.isfinite(global_index["S"]))

    def test_sensibilidad_registra_pasos_y_tolerancias(self):
        modelo = PIA.calcular_modelo(self.matriculas)
        sensibilidad = modelo["sensibilidad"]
        self.assertEqual(
            [dato["h"] for dato in sensibilidad["pasos"]],
            [0.1, 0.05, 0.025],
        )
        self.assertEqual(
            [dato["tolerancia"] for dato in sensibilidad["tolerancias"]],
            [1e-6, 1e-8, 1e-10],
        )
        self.assertTrue(all(math.isfinite(dato["S"]) for dato in sensibilidad["pasos"]))

    def test_informe_contiene_resultados_graficas_y_supuestos(self):
        modelo = PIA.calcular_modelo(self.matriculas)
        informe = PIA.crear_reporte_html(modelo, self.matriculas)
        self.assertIn("<svg", informe)
        self.assertIn("Factor de normalización térmica", informe)
        self.assertIn("En palabras simples", informe)
        self.assertIn("Qué significa", informe)
        self.assertIn("Suma ponderada (S base)", informe)
        self.assertIn("Documentos B–G — Resultados explicados", informe)
        self.assertIn("theta(0) = 25", informe)
        self.assertIn("Sensibilidad a la tolerancia", informe)
        self.assertIn("Código fuente completo", informe)
        self.assertIn("class Aplicacion", informe)
        self.assertIn("Validaciones numéricas", informe)


if __name__ == "__main__":
    unittest.main()
