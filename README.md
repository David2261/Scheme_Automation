# Scheme_Automation

## Description
*Examples of automation functional diagram, block diagram, tools specifications, stability factor, quality factor, switchboard layout, operator's workplace arrangement and presentation.*


## Tools

* [Kompass (КОМПАС-3D v20)](https://kompas.ru/)
* Program

*All calculations and all text are written in Russian, by [me](https://github.com/David2261 "Bulat Nasyrov")*

## Standards

| Purpose | Number |
| --- | ----------- |
| Functional scheme of automation | ГОСТ 21.404-85 |
| External connection diagram | СНиП 3.05.05-84 |


## Example
```pascal
IF "Start" THEN
	"RawFlow" := TRUE;
ELSE
	"RawFlow" := FALSE;
END_IF;

// Цепь 2
IF "Level_1" THEN
	"RawFlow" := TRUE;
	"OnOven" := FALSE;
ELSE
	"OnOven" := TRUE;
	"RawFlow" := FALSE;
END_IF;
```