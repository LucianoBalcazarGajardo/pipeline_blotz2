#!/usr/bin/env python3
import os
import glob
import json
import warnings
import argparse
import pandas as pd
import matplotlib.pyplot as plt

from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio.PDB import MMCIFParser, NeighborSearch

DISTANCIA_MINIMA = 2.0
DISTANCIA_MAXIMA = 3.5

RESIDUOS_PROTEINA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
RESIDUOS_ARN = {'A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC', 'DA', 'DT', 'DG', 'DC'}

def calcular_interacciones_cif(ruta_cif, dist_min, dist_max):
    parser = MMCIFParser(QUIET=True)
    estructura = parser.get_structure('modelo', ruta_cif)

    atomos_proteina, atomos_arn = [], []
    for modelo in estructura:
        for cadena in modelo:
            for residuo in cadena:
                nombre_res = residuo.get_resname().strip().upper()
                if nombre_res in RESIDUOS_PROTEINA:
                    atomos_proteina.extend(residuo.get_atoms())
                elif nombre_res in RESIDUOS_ARN:
                    atomos_arn.extend(residuo.get_atoms())

    if not atomos_proteina or not atomos_arn:
        return 0, 0, {}

    motor_busqueda = NeighborSearch(atomos_proteina)
    distancias_pares = {}
    choques_estericos = set()

    for atomo_arn in atomos_arn:
        atomos_cercanos = motor_busqueda.search(atomo_arn.get_coord(), dist_max)
        for atomo_prot in atomos_cercanos:
            distancia_real = atomo_arn - atomo_prot
            res_prot, res_arn = atomo_prot.get_parent(), atomo_arn.get_parent()

            id_prot = f"{res_prot.get_resname().strip()} {res_prot.get_id()[1]}"
            id_arn = f"{res_arn.get_resname().strip()} {res_arn.get_id()[1]}"
            pareja = (id_prot, id_arn)

            if distancia_real < dist_min:
                choques_estericos.add(pareja)
            else:
                if pareja not in distancias_pares or distancia_real < distancias_pares[pareja]:
                    distancias_pares[pareja] = distancia_real

    interacciones_validas = {pareja: dist for pareja, dist in distancias_pares.items() if pareja not in choques_estericos}
    return len(interacciones_validas), len(choques_estericos), interacciones_validas

def ejecutar_pipeline(carpeta_raiz, umbral_iptm):
    print(f"=== INICIANDO PIPELINE ===")
    print(f"Analizando directorio: {carpeta_raiz}")
    print(f"Umbral ipTM configurado: {umbral_iptm}\n")

    print("1. Extrayendo datos de archivos JSON...")
    datos_indice = []
    archivos_json = glob.glob(os.path.join(carpeta_raiz, '**', '*.json'), recursive=True)

    for ruta_json in archivos_json:
        try:
            with open(ruta_json, 'r') as f:
                dj = json.load(f)

            ptm, iptm = dj.get('ptm'), dj.get('iptm')
            if ptm is None or iptm is None: continue

            dir_json = os.path.dirname(ruta_json)
            cifs = glob.glob(os.path.join(dir_json, '*.cif'))
            if not cifs: continue

            datos_indice.append({
                'Carpeta': os.path.basename(dir_json),
                'Nombre_CIF': os.path.basename(cifs[0]),
                'Ruta_CIF': cifs[0],
                'pTM': float(ptm),
                'ipTM': float(iptm)
            })
        except Exception as e:
            pass

    df_indice = pd.DataFrame(datos_indice)
    if df_indice.empty:
        print("No se encontraron modelos válidos.")
        return

    df_indice.to_csv(os.path.join(carpeta_raiz, "01_indice_general_modelos.csv"), index=False)

    print("2. Generando histogramas de confianza...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.hist(df_indice['pTM'], bins=40, color='#4C72B0', edgecolor='black', alpha=0.8)
    ax1.set_title('Distribución Global de pTM', weight='bold')

    ax2.hist(df_indice['ipTM'], bins=40, color='#55A868', edgecolor='black', alpha=0.8)
    ax2.axvline(umbral_iptm, color='red', linestyle='dashed', linewidth=2, label=f'Corte: {umbral_iptm}')
    ax2.set_title('Distribución Global de ipTM', weight='bold')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(carpeta_raiz, "02_histograma_inicial_ptm_iptm.png"), dpi=300)
    plt.close()

    print(f"3. Aplicando corte de ipTM >= {umbral_iptm}...")
    df_filtrado = df_indice[df_indice['ipTM'] >= umbral_iptm].copy()
    print(f"Modelos filtrados: {len(df_filtrado)}")

    if df_filtrado.empty: return

    df_filtrado.to_csv(os.path.join(carpeta_raiz, "03_indice_modelos_filtrados.csv"), index=False)

    plt.figure(figsize=(8, 5))
    plt.hist(df_filtrado['ipTM'], bins=20, color='#DD8452', edgecolor='black')
    plt.title(f'Distribución ipTM (Modelos Filtrados >= {umbral_iptm})', weight='bold')
    plt.savefig(os.path.join(carpeta_raiz, "04_histograma_filtrado_iptm.png"), dpi=300)
    plt.close()

    print("4. Analizando interacciones estructurales...")
    res_int, det_enl = [], []

    for idx, row in df_filtrado.iterrows():
        n_val, n_choq, dic_int = calcular_interacciones_cif(row['Ruta_CIF'], DISTANCIA_MINIMA, DISTANCIA_MAXIMA)

        res_int.append({**row.to_dict(), 'Interacciones_Validas': n_val, 'Choques_Descartados': n_choq})
        for (r_prot, r_arn), dist in dic_int.items():
            det_enl.append({'Carpeta': row['Carpeta'], 'Nombre_CIF': row['Nombre_CIF'],
                            'Residuo_Proteina': r_prot, 'Residuo_ARN': r_arn, 'Distancia_A': round(dist, 3)})

    df_resumen_int = pd.DataFrame(res_int)
    df_resumen_int.to_csv(os.path.join(carpeta_raiz, "05_resumen_interacciones.csv"), index=False)

    if det_enl:
        pd.DataFrame(det_enl).to_csv(os.path.join(carpeta_raiz, "06_detalles_distancias.csv"), index=False)

    print("5. Generando histograma de interacciones...")
    plt.figure(figsize=(10, 6))
    max_int = df_resumen_int['Interacciones_Validas'].max()

    plt.hist(df_resumen_int['Interacciones_Validas'], bins=range(0, (max_int or 10) + 2), color='#8172B3', edgecolor='black', align='left')
    plt.title(f'Interacciones Válidas por Modelo ({DISTANCIA_MINIMA}Å - {DISTANCIA_MAXIMA}Å)', weight='bold')
    plt.savefig(os.path.join(carpeta_raiz, "07_histograma_interacciones.png"), dpi=300)
    plt.close()

    print("\n=== PIPELINE FINALIZADO CORRECTAMENTE ===")
    print(f"Los resultados se han guardado dentro de: {os.path.abspath(carpeta_raiz)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline de Análisis Blitz-2 para complejos Proteína-ARN")

    parser.add_argument("-i", "--input", required=True, help="Ruta a la carpeta que contiene las simulaciones a analizar.")
    parser.add_argument("-u", "--umbral", type=float, default=0.8, help="Valor de corte para ipTM (por defecto: 0.8)")

    # Volvemos a leer desde la consola
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"Error: El directorio '{args.input}' no existe.")
    else:
        ejecutar_pipeline(args.input, args.umbral)
