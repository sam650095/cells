async function processbtn(event) {
    event.preventDefault();
    const csrftoken = getCookie("csrftoken");
    
    toggleLoading(true, "processbutton");
    document.getElementById("imgbox").classList.remove("hidden");
    document.getElementById("nextbtn").classList.remove("hidden");
    const formData = new FormData();
    formData.append("file", selectedFile);
  
    const spatialanalysis_result = await fetchAPI(
      "/api/spatial",
      formData,
      csrftoken
    );
    document.getElementById("spatialanalysis_result").textContent =
      spatialanalysis_result.phenotyping_result;
    toggleLoading(false, "processbutton");
    
}