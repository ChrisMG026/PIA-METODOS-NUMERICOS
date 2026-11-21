import unittest

import PIA


class PruebasMetodosNumericos(unittest.TestCase):
    def setUp(self):
        identificaciones = ("2022022", "2105057", "2095958", "1950199")
        self.matriculas = {
            f"M{indice}": [int(digito) for digito in texto]
            for indice, texto in enumerate(sorted(identificaciones), start=1)
        }

    def test_interpolacion_elige_formula_compatible_con_la_malla(self):
        valor_forward = PIA.interpolacion_newton_adelante(
            4.5, (2, 4, 6, 8), (10.5, 18.7, 31.1, 47.3)
        )
        self.assertAlmostEqual(valor_forward, 21.421875, places=10)
        self.assertAlmostEqual(
            PIA.interpolacion_lagrange(2.5, (1, 2, 3, 4), (2.8, 4.6, 7.9, 12.5)),
            PIA.interpolacion_newton_adelante(
                2.5, (1, 2, 3, 4), (2.8, 4.6, 7.9, 12.5)
            ),
            places=10,
        )
        with self.assertRaises(ValueError):
            PIA.interpolacion_newton_adelante(2, (1, 2, 4), (1, 2, 3))

    def test_cs_usa_simpson_1_3_compuesto_y_expone_el_desglose(self):
        modelo = PIA.calcular_modelo(self.matriculas)
        c_s = modelo["documento_c"]["c_s"]
        self.assertAlmostEqual(c_s, 89.8 / 3, places=12)

        desglose = modelo["documento_c"]["desglose_c_s"]
        self.assertEqual(desglose["nodos"], (0, 1, 2, 3, 4))
        self.assertEqual(desglose["valores"], (3.0, 4.4, 6.9, 10.1, 15.0))
        self.assertEqual(desglose["num_intervalos"], 4)
        self.assertEqual(desglose["h"], 1)
        self.assertAlmostEqual(desglose["suma_extremos"], 18.0, places=12)
        self.assertAlmostEqual(desglose["suma_impares"], 14.5, places=12)
        self.assertAlmostEqual(desglose["suma_pares"], 6.9, places=12)
        self.assertAlmostEqual(desglose["suma_ponderada"], 89.8, places=12)
        self.assertAlmostEqual(desglose["resultado"], c_s, places=12)
        self.assertAlmostEqual(
            PIA.simpson_1_3_tabla(desglose["valores"], desglose["h"]),
            c_s,
            places=12,
        )

    def test_desglose_de_cs_aparece_en_html_y_en_operaciones(self):
        modelo = PIA.calcular_modelo(self.matriculas)
        informe = PIA.crear_reporte_html(modelo, self.matriculas)
        for texto in (
            "desglose de c_s con Simpson 1/3 compuesto",
            "Σ extremos (l0 + l4)",
            "89.8",
            "29.93333333",
            "Simpson 1/3 compuesto (sobre funciones y sobre tablas",
        ):
            self.assertIn(texto, informe)

        import tkinter as tk

        try:
            aplicacion = PIA.Aplicacion()
        except tk.TclError:
            self.skipTest("Entorno sin interfaz gráfica disponible")
        try:
            aplicacion.matriculas = self.matriculas
            aplicacion.modelo = modelo
            aplicacion.mostrar_operaciones()
            area = None
            for widget in aplicacion.content.winfo_children():
                for hijo in widget.winfo_children():
                    if isinstance(hijo, tk.Text):
                        area = hijo
            self.assertIsNotNone(area)
            contenido = area.get("1.0", "end")
            for texto in (
                "c_s: Simpson 1/3 compuesto",
                "Σ ponderada = 18 + 4·14.5 + 2·6.9 = 89.8",
                "c_s = (h/3)·Σ ponderada = (1/3)·89.8 = 29.93333333",
            ):
                self.assertIn(texto, contenido)
        finally:
            aplicacion.destroy()

    def test_ar_usa_la_menor_raiz_positiva_de_su_cubica(self):
        datos_b = PIA.calcular_coeficientes_documento_b()
        raices_positivas = [raiz for raiz in datos_b["raices_ar"] if raiz > 0]

        self.assertEqual(len(raices_positivas), 3)
        self.assertAlmostEqual(datos_b["Ar"], min(raices_positivas), places=12)

    def test_d_e_y_f_tienen_dependencia_y_metodos_declarados(self):
        modelo = PIA.calcular_modelo(self.matriculas)
        datos_d = modelo["documento_d"]
        datos_e = modelo["documento_e"]
        datos_f = modelo["documento_f"]

        self.assertAlmostEqual(
            datos_d["L3"],
            datos_d["k"] * 3**4 + datos_d["m"] * 3**3
            + datos_d["n"] * 3**2 + datos_d["p"] * 3 + datos_d["q"],
        )
        self.assertEqual(datos_e["kappa"], datos_d["L3"])
        self.assertAlmostEqual(datos_e["omega"], 12.0, places=12)
        self.assertEqual(datos_f["theta0"], 25.0)
        self.assertNotIn("escala_fuente", datos_f)
        self.assertGreater(datos_f["fuente_minima"], datos_f["capacidad_maxima"])
        self.assertEqual(datos_f["estado"], "divergente")
        self.assertIsNone(datos_f["theta_final"])
        self.assertLess(datos_f["tiempo_fallo"], 3.0)
        self.assertIn("Y/xi se mantiene por encima", datos_f["causa"])
        self.assertIn("capacidad máxima del término logístico", datos_f["causa"])

    def test_g_no_inventa_indice_si_f_no_es_finito(self):
        modelo = PIA.calcular_modelo(self.matriculas)
        datos_g = modelo["documento_g"]

        self.assertEqual(datos_g["estado"], "no disponible")
        self.assertNotIn("S", datos_g)
        self.assertIn("No es válido calcular", datos_g["causa"])

    def test_sensibilidad_registra_fallos_y_no_disponibilidad(self):
        sensibilidad = PIA.calcular_modelo(self.matriculas)["sensibilidad"]

        self.assertEqual(
            [dato["h"] for dato in sensibilidad["pasos"]],
            [0.1, 0.05, 0.025],
        )
        self.assertEqual(
            [dato["tolerancia"] for dato in sensibilidad["tolerancias"]],
            [1e-6, 1e-8, 1e-10],
        )
        self.assertTrue(all(dato["S"] is None for dato in sensibilidad["pasos"]))
        self.assertTrue(all(dato["theta_final"] is None for dato in sensibilidad["pasos"]))
        self.assertTrue(all(dato["S"] is None for dato in sensibilidad["tolerancias"]))

    def test_reporte_muestra_modelo_original_y_estado_parcial(self):
        modelo = PIA.calcular_modelo(self.matriculas)
        informe = PIA.crear_reporte_html(modelo, self.matriculas)

        self.assertIn("<svg", informe)
        self.assertIn("Interpretación", informe)
        self.assertIn("Qué significa", informe)
        self.assertIn("Documentos B–G · parámetros reconstruidos", informe)
        self.assertIn("theta(0)=25", informe)
        self.assertIn("Divergente", informe)
        self.assertIn("capacidad máxima del término logístico", informe)
        self.assertIn("sin normalización", informe)
        self.assertIn("No disponible", informe)
        self.assertIn("Runge-Kutta de cuarto orden", informe)
        self.assertIn('class="grid"', informe)
        self.assertIn('class="badge badge-danger"', informe)
        self.assertIn("Escala recortada", informe)
        self.assertIn("page-break-inside:avoid", informe)
        self.assertNotIn("Factor de normalización térmica", informe)
        self.assertNotIn("omega = kappa", informe)
        self.assertIn("Código fuente completo", informe)
        self.assertIn("Interpretación y validación", informe)


if __name__ == "__main__":
    unittest.main()
