// Near real-time feed for poll/index.html (no "Publish to web" cache delay).
// 1. In the Google Sheet linked to the form: Extensions > Apps Script, paste this file, Save.
// 2. Deploy > New deployment > Type: Web app > Execute as: Me > Who has access: Anyone > Deploy.
// 3. Copy the Web app URL (ends with /exec) into CSV_URL in poll/index.html.
function doGet() {
  const sheet = SpreadsheetApp.getActiveSpreadsheet().getSheets()[0];   // form responses sheet
  const rows = sheet.getDataRange().getDisplayValues();
  const csv = rows.map(r => r.map(c => '"' + String(c).replace(/"/g, '""') + '"').join(",")).join("\n");
  return ContentService.createTextOutput(csv).setMimeType(ContentService.MimeType.CSV);
}
