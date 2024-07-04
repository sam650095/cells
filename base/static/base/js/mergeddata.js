// proccess button click
async function processbtn(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  const process_result = await fetchAPI("/api/merge", 0, csrftoken);
  console.log(process_result);
  // addata_results div
  let adata_results = document.getElementById("adata_results");
  adata_results.innerHTML = "";
  let ul = document.createElement("ul");
  ul.classList.add("space-y-1", "text-gray-500", "list-disc", "list-inside");
  process_result["adata_results"].forEach((item) => {
    let li = document.createElement("li");
    li.textContent = item;
    ul.appendChild(li);
  });
  adata_results.appendChild(ul);
  // next page btn
  document.getElementById("nextbtn").classList.remove("hidden");
}
// fetch api
async function fetchAPI(url, formData, csrftoken) {
  const response = await fetch(url, {
    method: "POST",
    body: formData,
    headers: {
      "X-CSRFToken": csrftoken,
    },
  });

  if (!response.ok) {
    throw new Error("Network response was not ok");
  }

  return response.json();
}