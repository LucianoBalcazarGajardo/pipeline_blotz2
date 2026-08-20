#!/usr/bin/env python3
import os
import re
import pandas as pd
import subprocess
import argparse
from Bio import SeqIO

# ==========================================
# Ruta fija al indice con fasta de todos los ARN en estudio 
# Editar esta ruta en el ambiente de trabajo 
RUTA_FASTA_MAESTRO = "/ruta/absoluta/a/tu/Secuencias_ARN_combinaciones.fasta"
# ==========================================

def generar_weblogo_desde_fasta(carpeta, proteina, porcentaje):
    print(f"=== INICIANDO EXTRACCIÓN PARA WEBLOGO (Top {porcentaje}%) ===")
    
    # IMPORTANTE: Con el archivo 05 con el nombre de la proteína incluido
    ruta_csv = os.path.join(carpeta, f"05_resumen_interacciones_{proteina}.csv")
    
    if not os.path.exists(ruta_csv):
        print(f"Error: No se encontró {ruta_csv}. Ejecuta el pipeline de filtrado primero con la misma proteína.")
        return

    # 1. Leer y aislar los mejores modelos
    df = pd.read_csv(ruta_csv)
    df_ordenado = df.sort_values(by='Interacciones_Validas', ascending=False).reset_index(drop=True)
    cantidad_top = max(1, int(len(df_ordenado) * (porcentaje / 100.0)))
    df_top = df_ordenado.head(cantidad_top)
    
    print(f"1. Aislados los {cantidad_top} mejores modelos del CSV.")

    # 2. Cargar el FASTA maestro
    diccionario_arn = {}
    try:
        for record in SeqIO.parse(RUTA_FASTA_MAESTRO, "fasta"):
            diccionario_arn[record.id] = str(record.seq)
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo FASTA maestro en {RUTA_FASTA_MAESTRO}")
        print("Verificar que la variable RUTA_FASTA_MAESTRO sea correcta.")
        return
        
    nombre_base = f"top{porcentaje}pct_{proteina}"
    ruta_csv_out = os.path.join(carpeta, f"{nombre_base}_mejores_modelos.csv")
    ruta_fasta_out = os.path.join(carpeta, f"{nombre_base}_secuencias.fasta")
    ruta_logo = os.path.join(carpeta, f"{nombre_base}_weblogo.png")
    
    df_top.to_csv(ruta_csv_out, index=False)
    
    # 3. Cruzar datos
    print("2. Cruzando identificadores y generando nuevo FASTA...")
    encontrados = 0
    with open(ruta_fasta_out, 'w') as f_out:
        for idx, row in df_top.iterrows():
            # Busca "ARN" seguido de números en la ruta
            match = re.search(r'(ARN\d+)', str(row['Carpeta']) + str(row['Nombre_CIF']))
            if match:
                arn_id = match.group(1)
                if arn_id in diccionario_arn:
                    f_out.write(f">{arn_id}\n{diccionario_arn[arn_id]}\n")
                    encontrados += 1

    if encontrados == 0:
        print(" No se emparejó ningún nombre del CSV con los IDs del FASTA maestro.")
        return

    # 4. Generar el gráfico
    print("3. Ejecutando WebLogo3 headless")
    try:
        comando = [
            "weblogo", 
            "-f", ruta_fasta_out, 
            "-o", ruta_logo, 
            "-F", "png", 
            "--resolution", "300", 
            "--title", f"Top {porcentaje}% - {proteina}",  # Título más corto
            "--title-fontsize", "10",                    # NUEVO: Achica la letra del título (por defecto es 12)
            "--size", "large",                           # NUEVO: Agranda las proporciones de toda la imagen
            "--units", "probability",    
            "--sequence-type", "rna"     
        ]
        subprocess.run(comando, check=True)
        subprocess.run(comando, check=True)
        print(f"=== PROCESO TERMINADO ===")
        print(f"Archivos guardados:")
        print(f" - {ruta_csv_out}")
        print(f" - {ruta_fasta_out}")
        print(f" - {ruta_logo}")
    except FileNotFoundError:
        print("Error: WebLogo no está instalado. Ejecuta 'pip install weblogo' en tu ambiente virtual.")
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar WebLogo: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Análisis secundario: WebLogo de los mejores ARNs")
    parser.add_argument("-i", "--input", required=True, help="Ruta a la carpeta de resultados")
    parser.add_argument("-p", "--proteina", required=True, help="Nombre de la proteína")
    parser.add_argument("-pct", "--porcentaje", type=float, required=True, help="Porcentaje de mejores modelos (ej. 10 para 10%)")
    
    args = parser.parse_args()
    
    generar_weblogo_desde_fasta(args.input, args.proteina, args.porcentaje)
