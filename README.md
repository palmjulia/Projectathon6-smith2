# Selectanfrage für den 6. Projectathon der MII: SMITH
Datum: 25.11.21

Autorin: [julia.palm@med.uni-jena.de](mailto:julia.palm@med.uni-jena.de).

Dieses Project führt die Select-Anfrage für das SMITH Projekt im Rahmen des 6. Projectathons aus. Hier ist eine zentrale Analyse vorgesehen. Dafür erzeugt dieses Skript zwei Tabellen mit den für die Analyse benötigten Inhalten. Diese Tabellen sollen zentral zusammengeführt und an die datenauswertendende Stelle übergeben werden.

Das Readme beschreibt zunächst die technischen Details der Verwendung. Darunter sind die verwendeten CodeSysteme/Ressourcen/Profile und der konzeptionelle Ablauf der Abfrage beschrieben.

## Verwendung
Es gibt zwei Möglichkeiten diese R-Skripte auszuführen: Direkt in R oder in einem Docker Container. Beide werden im folgenden beschrieben.

### Ausführung in R
#### Vor der ersten Nutzung
1. Um die Selectanfrage durchzuführen, muss der Inhalt des Git-Repository auf einen Rechner (PC, Server) gezogen werden, von dem aus der REST-Endpunkt des gewünschten FHIR-Servers (z.B. FHIR-Server der Clinical Domain im DIZ) erreichbar ist. 

2. Auf diesem Rechner muss R (aber nicht notwendigerweise RStudio) als genutzte Laufzeitumgebung installiert sein.

3. Die mitgelieferte Datei `./config.R.default` muss nach `./config.R` kopiert werden und lokal angepasst werden (FHIR-Endpunkt, ggf. Authentifizierung, SSL peer verification); Erklärungen dazu finden sich direkt in dieser Datei. Eine Authentifizierung mit Basic Authentication oder Bearer Token ist möglich. Dafür müssen in `config.R` die Variable `authentication` und die zugehörigen Zugangsdaten (`password`/`username` bzw. `token`) angepasst werden.
Wenn die Abfrage auf einem Server laufen sollen, der sowohl konsentierte als auch nicht konsentierte Daten enthält, so kann durch setzen der Variable `filterConsent <- TRUE` dafür gesorgt werden, dass nur Daten von Patienten extrahiert werden, auf die eine Consent-Ressource mit einem `2.16.840.1.113883.3.1937.777.24.5.3.8` (*MDAT_wissenschaftlich_nutzen_EU_DSGVO_NIVEAU*) Code verweist.  

4. Wenn die App über `runSmith_select.bat` (unter Windows) gestartet soll, muss in dieser der Pfad zur Datei `Rscript.exe` geprüft und ggf. angepasst werden (z.B. `C:\Program Files\R\R-4.0.4\bin\Rscript.exe`).


#### Start des Skripts
Beim ersten Start des Skripts wird überprüft, ob die zur Ausführung notwendigen R-Pakete (`fhircrackr`, `data.table`) vorhanden sind. Ist dies nicht der Fall, werden diese Pakete nachinstalliert – dieser Prozess kann einige Zeit in Anspruch nehmen.

##### Batch-Datei/Shell-Skript
**Unter Windows**: Mit der Batch-Datei `runSmith_select.bat`.
Beim ersten Ausführen sollte diese ggf. als Administrator gestartet werden (über Eingabeaufforderung oder Rechtsklick), wenn die ggf. notwendigen Berechtigungen zum Nachinstallieren der R-Pakete sonst nicht vorhanden sind. Nach der ersten Installation reicht dann ein Doppelklick zum Starten.

**Unter Linux**: Mit dem Shell-Skript `runSmith_selectr.sh`. Das Shell-Skript muss ausführbar sein und ggf. beim ersten Ausführen mittels `sudo` gestartet werden, wenn ein Nachinstallieren der R-Pakete außerhalb des User-Kontexts erforderlich ist.

#### R/RStudio
Durch Öffnen des R-Projektes (`Projectathon6-smith2.Rproj`) mit anschließendem Ausführen der Datei `smith_select.R` innerhalb von R/RStudio. Auch hier werden beim ersten Ausführen ggf. notwendige R-Pakete nachinstalliert.


## Ausführung im Docker Container
Um die Abfrage in einem Docker Container laufen zu lassen gibt es drei Möglichkeiten:

**A) Image von DockerHub ziehen:**
1. Git-Respository klonen: `git clone https://github.com/palmjulia/Projectathon6-smith2.git`
2. Verzeichniswechsel in das lokale Repository: `cd Projectathon6-smith2`
3. Konfiguration lokal anpassen: `./config.R.default` nach `./config.R` kopieren und anpassen 
4. Image downloaden und Container starten: `docker run --name projectathon6-smith2 -v "$(pwd)/errors:/errors" -v "$(pwd)/Bundles:/Bundles" -v "$(pwd)/Ergebnisse:/Ergebnisse" -v "$(pwd)/config.R:/config.R" palmjulia/projectathon6-smith2`


**B) Image bauen mit Docker Compose:**

1. Git-Respository klonen: `git clone https://github.com/palmjulia/Projectathon6-smith2.git`
2. Verzeichniswechsel in das lokale Repository: `cd Projectathon6-smith2`
3. Konfiguration lokal anpassen: `./config.R.default` nach `./config.R` kopieren und anpassen 
4. Image bauen und Container starten: `docker compose up -d`

Zum Stoppen des Containers `docker compose stop`. Um ihn erneut zu starten, `docker compose start`.

**C) Image bauen ohne Docker Compose**

1. Git-Respository klonen: `git clone https://github.com/palmjulia/Projectathon6-smith2.git`
2. Verzeichniswechsel in das lokale Repository: `cd Projectathon6-smith2`
3. Image bauen: `docker build -t projectathon6-smith2 .` 
4. Konfiguration lokal anpassen: `./config.R.default` nach `./config.R` kopieren und anpassen 
5. Container starten: `docker run --name projectathon6-smith2 -v "$(pwd)/errors:/errors" -v "$(pwd)/Bundles:/Bundles" -v "$(pwd)/Ergebnisse:/Ergebnisse" -v "$(pwd)/config.R:/config.R" projectathon6-smith2`

Erklärung:

-  `-v "$(pwd)/config.R:/config.R""` bindet die lokal veränderte Variante des config-Files ein. Wenn dieses geändert wird, reicht es, den Container neu zu starten (`docker stop Projectathon6-smith2`, config.R ändern, dann `docker start Projectathon6-smith2`), ein erneutes `docker build` ist nicht nötig.


-----------------------------------------------------------------------------------------------

## Selbstsignierte Server-Zertifikate

Falls der verwendete FHIR-Server ein selbst-signiertes Zertifikat für HTTPS verwendet, ist es notwendig das zugehörige Root-Zertifikat zu übergeben. R nutzt die Systemzertifikate des Betriebssystems hierfür.

Wenn die Abfrage direkt (ohne Docker) ausgeführt wird muss sich das Root-Zertifikat in den systemweit vertrauenswürdigen Zertifikaten befinden. Dies ist Betriebssystemabhängig. Für die Ausführung im Docker-Container kann das Root-Zertifikat als Volume eingebunden werden.

**Docker**
Angenommen das Zertifikat liegt unter `./localca.pem` muss folgender zusätzlicher Parameter beim `docker run ...`-Befehl übergeben werden:
`-v ./localca.pem:/usr/local/share/ca-certificates.crt`

**Docker-Compose**
Angenommen das Zertifikat liegt unter `./localca.pem` muss folgender Eintrag in der `docker-compose.yml` gemacht werden:
```
[...]
  volumes:
  - ./localca.pem:/usr/local/share/ca-certificates.crt
[...]
```

Hintergrund:
Das Image basiert auf Debian und führt bei jedem Start den Befehl `update-ca-certificates` aus, mit dem es die Zertifikate unter `/usr/local/share/ca-certificates` einliest und zu den Vertrauenswürdigen hinzufügt.


## Output
Das Skript erzeugt mehrere Ordner im Projekt-Directory. Um für den Projectathon eine möglichst einfache übersichtliche Lösung zu bekommen, werden alle files, die darin erzeugt werden bei mehrmaligem Ausführen ggf. einfach überschrieben.

### Ergebnisse
Wenn die Abfrage erfolgreich durchgeführt wurde, sind hier zwei semikolongetrennte csv-Dateien (= lässt sich durch Doppelklick in Excel öffnen) zu finden. Die enthaltenen Variablen sind hier kurz erklärt: 

**Kohorte.csv**

Diese Tabelle enthält eine Kombination von Informationen aus der Patient Ressource, der Encounter Ressource und der Observation Ressource. Sie enthält alle Fälle, für die es eine Observation mit einer NTproBNP Messung im geforderten Zeitraum gibt.

|Variable             | Bedeutung|
|---------------------|----------|
|subject              | Logical id der Patient Ressource|
|encounter.start      | Startzeitpunkt des Einrichtungskontakt-Encounters, der zeitlich zur NTproBNP-Messung gehört.|
|encounter.end        | Stoppzeitpunkt des Einrichtungskontakt-Encounters, der zeitlich zur NTproBNP-Messung gehört.|
|serviceType          | ServiceType (Stationsschlüssel) des Einrichtungskontakt-Encounters, der zeitlich zur NTproBNP-Messung gehört.|
|NTproBNP.date        | Datum (effectiveDateTime) der Observation mit der NTproBNP-Messung.|
|NTproBNP.value       | Wert (valueQuantity.value) der Observation mit der NTproBNP-Messung.|
|NTproBNP.unit        | Code der Einheit (valueQuantity.code) der NTproBNP-Messung.|
|NTproBNP.unitSystem  | Codesystem der Einheit (valueQuantity.code) der NTproBNP-Messung.|
|gender               | Geschlecht (gender) der Patient Ressource|
|birthdate            | Geburtsdatum (birthDate) der Patient Ressource|

**Diagnosen.csv**

Diese Tabelle enthält alle Diagnosen der Patienten aus Kohorte.csv, die einen der im Antrag genannten ICD-Kodes enthalten.

|Variable             | Bedeutung|
|---------------------|----------|
|clinicalStatus.code        |Code des ClinicalStatus|
|clinicalStatus.system      |Codesystem des ClinicalStatus|
|verificationStatus.code    |Code des verificationStatus|
|verificationStatus.system  |CodeSystem des verificationStatus|
|code                       |ICD-Code|
|code.system                |Codesystem des ICD Codes|
|subject                    |ID der zugehörigen Patient Ressource|
|onsetPeriod.start          |Onset Start|
|onsetPeriod.end            |Onset End|
|recordedDate               |recordedDate|

### Bundles
Dieser Ordner enthält die heruntergeladenen Bundles und kann der Kontrolle dienen, falls die Tabellen nicht aussehen wie erwartet.

### errors
Dieser Ordner enthält ggf. Fehlermeldungen, wenn die Abfrage nicht erfolgreich durchgeführt werden kann.

## Verwendete Codesysteme
Das Skript verwendet folgende Codesysteme:

- *http://loinc.org* für `Observation.code.coding.system` -> Dieses System wird für den Download per FHIR Search verwendet
- *xxxxxxicd-10xxxxx* für `Condition.code.coding.system` -> Im Skript wird per String matching nach irgendeinem CodeSystem gesucht, das die Zeichen *icd-10* enthält (Für Kompatibilität bei Änderungen des konkreten ICD-Systems)
- *urn:oid:2.16.840.1.113883.3.1937.777.24.5.3* für `Consent.provision.provision.code.coding.system` -> Nur verwendet, wenn konsentierte Daten über Consent Ressource selektiert werden, also wenn `filterConsent <- TRUE` in config.R.

## Verwendete Profile/Datenelemente
Die Abfragen sind auf Basis der MII Profile für die entsprechenden Ressouren geschrieben. Die Skripte sind kompatibel mit dem jeweils neuesten Release der verfügbaren Major Versionen (1 und, wenn vorhanden, 2). In den meisten Fällen werden die Skripte auch funktionieren, wenn nicht genau dieser Release des Profils verwendet wird. Es ist jedoch zwingend notwendig, dass die Ressourcen das jeweilige Profil im `Resource.meta.profile` Element benennen. Im folgenden ist für jeden verwendeten Ressourcentyp beschrieben, welche Elemente für die FHIR-Search-Abfrage an den Server verwendet werden (diese Elemente *müssen* vorhanden sein, damit kein Fehler geworfen wird) und welche Elemente im Skript extrahiert und in die Ergebnistabellen geschrieben werden. 

### Modul Labor: Observation
Profil: `https://www.medizininformatik-initiative.de/fhir/core/modul-labor/StructureDefinition/ObservationLab`

Version: [1.0.6](https://simplifier.net/packages/de.medizininformatikinitiative.kerndatensatz.laborbefund/1.0.6/~introduction)

Für Servabfrage verwendete Elemente:

- `Observation.subject`
- `Observation.code`
- `Observation.effective`
- `Observation.meta.profile`

Extrahierte Elemente:

- `Observation.effectiveDateTime`
- `Observation.subject.reference`
- `Observation.valueQuantity.value`
- `Observation.valueQuantity.code`
- `Observation.valueQuantity.system`

### Modul Person: Patient
Profil: `https://www.medizininformatik-initiative.de/fhir/core/modul-person/StructureDefinition/Patient`

Version: [2.0.0-alpha3](https://simplifier.net/packages/de.medizininformatikinitiative.kerndatensatz.person/2.0.0-alpha3) bzw. [1.0.14](https://simplifier.net/packages/de.medizininformatikinitiative.kerndatensatz.person/1.0.14)

Für Servabfrage verwendete Elemente:

- keine

Extrahierte Elemente:

- `Patient.id`
- `Patient.gender`
- `Patient.birthDate`

### Modul Fall: Encounter
Profil: `https://www.medizininformatik-initiative.de/fhir/core/modul-fall/StructureDefinition/KontaktGesundheitseinrichtung`

Version: [1.0.1](https://simplifier.net/packages/de.medizininformatikinitiative.kerndatensatz.fall/1.0.1)

Für Servabfrage verwendete Elemente:

- `Encounter.meta.profile`
- `Encounter.subject.reference`

Extrahierte Elemente:

- `Encounter.subject.reference`
- `Encounter.period.start `
- `Encounter.period.end`
- `Encounter.serviceType`

### Modul Diagnose: Condition
Profil: `https://www.medizininformatik-initiative.de/fhir/core/modul-diagnose/StructureDefinition/Diagnose`

Version: [2.0.0-alpha3](https://simplifier.net/packages/de.medizininformatikinitiative.kerndatensatz.diagnose/2.0.0-alpha3) bzw. [1.0.4](https://simplifier.net/packages/de.medizininformatikinitiative.kerndatensatz.diagnose/1.0.4)

Für Servabfrage verwendete Elemente:

- `Condition.meta.profile`
- `Condition.subject.reference`

Extrahierte Elemente:

- `Condition.clinicalStatus.coding.code`
- `Condition.clinicalStatus.coding.system`
- `Condition.verificationStatus.coding.code`
- `Condition.verificationStatus.coding.system`
- `Condition.code.coding.code`
- `Condition.code.coding.system`
- `Condition.subject.reference`
- `Condition.onsetPeriod.start`
- `Condition.onsetPeriod.end`
- `Condition.recordedDate`

### Modul Consent: Consent
Wird nur verwendet, wenn `filterConsent <- TRUE` in `config.R`, d.h. wenn konsentierte und nicht konsentierte Daten durch das R-Skript gefiltert werden müssen und die Trennung nicht durch unterschiedliche FHIR-Repositories erfolgt.

Profil: `http://fhir.de/ConsentManagement/StructureDefinition/Consent`
Version: [1.0.0](https://simplifier.net/packages/de.einwilligungsmanagement/1.0.0)
Beispiel: https://simplifier.net/packages/de.einwilligungsmanagement/1.0.0/files/405292

Für Servabfrage verwendete Elemente:

- `Consent.patient.reference`

Extrahierte Elemente:

- `Consent.patient.reference`
- `Consent.provision.provision.code.coding.code`
- `Consent.provision.provision.code.coding.system`

## Konzeptioneller Ablauf der Abfrage
Prinzipiell geht das Skript wie folgt vor:

1) Lade alle Observations, die eine NTproBNP-Messung darstellen (`loinc code` ist einer aus: `33763-4,71425-3,33762-6,83107-3, 83108-1`), welche ein `effective` im Zeitraum `2019-01-01 - 2021-12-31` haben und das Profil `https://www.medizininformatik-initiative.de/fhir/core/modul-labor/StructureDefinition/ObservationLab` erfüllen, herunter. Ziehe dazu außerdem (mit `_include`) die zugehörigen Patient Ressourcen.

Request: `"[base]/Observation?code=http://loinc.org%7C33763-4,http://loinc.org%7C71425-3,http://loinc.org%7C33762-6,http://loinc.org%7C83107-3,http://loinc.org%7C83108-1&date=ge2019-01-01&date=le2021-12-31&_include=Observation:patient&_profile=https://www.medizininformatik-initiative.de/fhir/core/modul-labor/StructureDefinition/ObservationLab"`

2) Extrahiere die IDs der Patient Ressourcen und verwende sie um alle zugehörigen Encounter herunterzuladen. Dafür werden die Abfragen soweit gesplittet, dass der jeweilige GET-Request eine Länge von 1800 Zeichen nicht überschreitet.(Mir ist klar, das POST eleganter ist, aber wir versuchen die Anforderungen an den Server/Rechte möglichst niedrig zu halten.) Die Encounter müssen das Profil KontaktGesundheitseinrichtung der MII implementieren. 

Request (beispielhaft für Patient ids `xxx` und `yyy`): `[base]/Encounter?_profile=https%3A%2F%2Fwww.medizininformatik-initiative.de%2Ffhir%2Fcore%2Fmodul-fall%2FStructureDefinition%2FKontaktGesundheitseinrichtung&subject=xxx%2Cyyy`

3) Filtere die Encounter, sodass nur die Encounter übrig bleiben, in deren `period` eine NTproBNP Observation liegt.

4) Lade nach dem gleichen Prinzip wie in 2) alle Condition Ressourcen herunter, die auf die in 1) selektierten Patienten verweisen. Die Conditions müssen das Profil Diagnose der MII implementieren. 

Request (beispielhaft für Patient ids `xxx` und `yyy`):
`[base]/Condition?_profile=https%3A%2F%2Fwww.medizininformatik-initiative.de%2Ffhir%2Fcore%2Fmodul-diagnose%2FStructureDefinition%2FDiagnose&subject=xxx%2Cyyy`

5) Filtere die Diagnosen so, dass nur noch Diagnosen enthalten sind, die einen der im Antrag genannten ICD-Codes enthalten.


