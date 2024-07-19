const fileInput = document.getElementById("file_upload");
const uploadedDiv = document.getElementById("uploaded");
let selectedFile = null;
fileInput.addEventListener("change", function (event) {
  const file = event.target.files[0];
  if (file) {
    uploadedDiv.textContent = file.name;
    selectedFile = file;
  } else {
    uploadedDiv.textContent = "";
    selectedFile = null;
  }
});
async function processbtn(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  toggleLoading(true, "processbutton");
  const formData = new FormData();
  formData.append("file", selectedFile);

  const identifythegates_result = await fetchAPI(
    "/api/phenotyping",
    formData,
    csrftoken
  );
  console.log(identifythegates_result);
  toggleLoading(false, "processbutton");

  //   document.getElementById("btnbox").classList.remove("hidden");
  //   document.getElementById("nextbtn").classList.remove("hidden");
}
